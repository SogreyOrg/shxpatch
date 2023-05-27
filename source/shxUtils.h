#pragma once
#include<string>
#include <iostream>
#include <memory>
using namespace std;

class shxUtils
{
public:


    template< typename... Args >
    static std::string string_format(const char* format, Args... args)
    {
        size_t length = std::snprintf(nullptr, 0, format, args...);
        if (length <= 0)
        {
            return "";
        }

        char* buf = new char[length + 1];
        std::snprintf(buf, length + 1, format, args...);

        std::string str(buf);
        delete[] buf;
        return std::move(str);
    }
};