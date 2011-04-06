#ifndef APP_INCLUDED_H
#define APP_INCLUDED_H 1

#include "common/app-base.hpp"

class App : public AppBase {

public:

    App(int debug,
        int profile,
        int verbose,
        std::vector<int> const & device_list)
        : AppBase(debug, profile, verbose),
          device_list_m(device_list)
        {
        }
    
    virtual void host_run();
	
private:

    virtual std::string const & get_device_program_text();

    std::vector<int> const & device_list_m;

};

#endif
