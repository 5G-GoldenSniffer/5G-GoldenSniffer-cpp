#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <pthread.h>
#include "util.h"

bool verbose = false;

void save_fc32_debug(const cfloat* data, size_t length, const std::string& base_name)
{
    static int counter = 0;
    std::string filename;
    do {
        std::ostringstream oss;
        oss << base_name << "-" << std::setw(4) << std::setfill('0') << counter << ".fc32";
        filename = oss.str();
        counter++;
    } while (file_exists(filename));
    std::ofstream file(filename, std::ios::binary);
    if (file.is_open()) {
        file.write(reinterpret_cast<const char*>(data), length * sizeof(cfloat));
        file.close();
    }
}

void save_fc32_debug(const std::vector<cfloat>& data, const std::string& base_name)
{
    save_fc32_debug(data.data(), data.size(), base_name);
}

cpu_set_t rx_cpus;
cpu_set_t worker_cpus;

bool parse_cpu_list(const std::string& list, cpu_set_t& cpus, std::string& err)
{
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if(sched_getaffinity(0, sizeof(allowed), &allowed) != 0){
        err = "sched_getaffinity failed";
        return false;
    }
    CPU_ZERO(&cpus);
    if(list == "all"){
        cpus = allowed;
        return true;
    }
    std::vector<int> dropped;
    std::istringstream ss(list);
    std::string tok;
    while(std::getline(ss, tok, ',')){
        // a, a-b or a-b:stride
        int a, b, stride = 1;
        char c1 = 0, c2 = 0, extra;
        int n = sscanf(tok.c_str(), " %d %c %d %c %d %c", &a, &c1, &b, &c2, &stride, &extra);
        if(n == 1)
            b = a;
        else if(!((n == 3 && c1 == '-') || (n == 5 && c1 == '-' && c2 == ':'))){
            err = "\"" + tok + "\" is not a CPU, a range a-b or a range with a stride a-b:s";
            return false;
        }
        if(a < 0 || b < a || stride < 1 || b >= CPU_SETSIZE){
            err = "\"" + tok + "\" is not a valid CPU range";
            return false;
        }
        for(int i = a; i <= b; i += stride){
            if(CPU_ISSET(i, &allowed))
                CPU_SET(i, &cpus);
            else
                dropped.push_back(i);
        }
    }
    if(!dropped.empty()){
        std::cout << "Warning: CPU list \"" << list << "\": CPUs not available to this process, ignored:";
        for(int i : dropped)
            std::cout << " " << i;
        std::cout << std::endl;
    }
    if(CPU_COUNT(&cpus) == 0){
        err = "no usable CPU in \"" + list + "\"";
        return false;
    }
    return true;
}

std::string cpu_list_string(const cpu_set_t& cpus)
{
    std::ostringstream oss;
    for(int i = 0; i < CPU_SETSIZE; i++)
        if(CPU_ISSET(i, &cpus))
            oss << (oss.tellp() ? "," : "") << i;
    return oss.str();
}

void pin_self(const cpu_set_t& cpus)
{
    int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpus), &cpus);
    if(rc != 0)
        std::cout << "Warning: pthread_setaffinity_np failed (" << rc << ")" << std::endl;
}
