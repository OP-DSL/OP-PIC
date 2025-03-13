#pragma once

#include "opp_defs.h"

namespace opp
{
    class Params
    {
        public:
            Params(std::string fileanme);

            template <typename T>
            inline void add(std::string key, T value) // refer specilized template functions
            {
                UNUSED_VAR(value);
                std::cerr << "Params add not implemented for type of : " << key << std::endl;
                // exit(-1);
            }

            template <typename T>
            inline T get(std::string key, bool error = true) // refer specilized template functions
            {
                std::cerr << "Params get not implemented for type of : " << key << std::endl;
                // exit(-1);
                return NULL;
            }

            void write(std::ostream &out);

        private:
            std::map<std::string, OPP_STRING> str_params;
            std::map<std::string, OPP_INT> int_params;
            std::map<std::string, OPP_REAL> real_params;
            std::map<std::string, OPP_BOOL> bool_params;
    };
    
    template <>
    inline void Params::add<OPP_INT>(std::string key, OPP_INT value)
    {
        int_params[key] = value;
    }

    template <>
    inline void Params::add<OPP_REAL>(std::string key, OPP_REAL value)
    {
        real_params[key] = value;
    }

    template <>
    inline void Params::add<OPP_BOOL>(std::string key, OPP_BOOL value)
    {
        bool_params[key] = value;
    }

    template <>
    inline void Params::add<OPP_STRING>(std::string key, OPP_STRING value)
    {
        str_params[key] = value;
    }

    template <>
    inline OPP_INT Params::get<OPP_INT>(std::string key, bool error)
    {
        auto a = int_params.find(key);
        if (a != int_params.end()) return a->second;
        if (error)
            std::cerr << "ERROR : " << key << " (OPP_INT) not found in configs" << std::endl;
        return INT_MAX; 
    }

    template <>
    inline OPP_REAL Params::get<OPP_REAL>(std::string key, bool error)
    {
        auto a = real_params.find(key);
        if (a != real_params.end()) return a->second;
        if (error)
            std::cerr << "ERROR : " << key << " (OPP_REAL) not found in configs" << std::endl;
        return DBL_MAX; 
    }
    
    template <>
    inline OPP_BOOL Params::get<OPP_BOOL>(std::string key, bool error)
    {
        auto a = bool_params.find(key);
        if (a != bool_params.end()) return a->second;
        if (error)
            std::cerr << "ERROR : " << key << " (BOOL) not found in configs" << std::endl;
        return false; 
    }

    template <>
    inline OPP_STRING Params::get<OPP_STRING>(std::string key, bool error)
    {
        auto a = str_params.find(key);
        if (a != str_params.end()) return a->second;
        if (error)
            std::cerr << "ERROR : " << key << " (OPP_STRING) not found in configs" << std::endl;
        return "NOT FOUND"; 
    }
};
