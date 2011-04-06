#ifndef APP_BASE_INCLUDED_H
#define APP_BASE_INCLUDED_H 1

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cstdio>
#include <cstdlib>

#define __CL_ENABLE_EXCEPTIONS
#include "common/cl.hpp"

///
// A class providing OpenCL boiler plate for simple applications
///
class AppBase {

public:

    AppBase(int debug,
            int profile,
            int verbose)
        : debug_m(debug),
          profile_m(profile),
          verbose_m(verbose)
        {
            int rc = CL_SUCCESS;

            // Check we have an OpenCL implementation
            cl::Platform::get(&platform_m);
            if (platform_m.size() == 0) {
                std::cerr << "Platform size 0\n";
                throw;
            }
            cl_context_properties context_properties[] = {
                CL_CONTEXT_PLATFORM,
                (cl_context_properties)(platform_m[0])(),
                0
            };
            context_m = cl::Context(CL_DEVICE_TYPE_ALL, context_properties);

            // Get devices, and create command queues for them
            device_m = context_m.getInfo<CL_CONTEXT_DEVICES>();
            if (verbose_m) {
                std::cerr << "[System]\n";
                std::cerr << "  number of devices: " << device_m.size() << "\n";
                for (size_t i = 0; i < device_m.size(); ++i) {
                    cl_device_type type = device_m[i].getInfo<CL_DEVICE_TYPE>();
                    std::string device_type;
                    if (type == CL_DEVICE_TYPE_CPU) {
                        device_type = std::string("CPU");
                    } else if (type == CL_DEVICE_TYPE_GPU) {
                        device_type = std::string("GPU");
                    } else if (type == CL_DEVICE_TYPE_ACCELERATOR) {
                        device_type = std::string("ACCELERATOR");
                    } else {
                        device_type = std::string("Huh?");
                    }
                    cl_uint max_work_item_dimensions = device_m[i].getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>();
                    std::ostringstream max_work_item_sizes;
                    for (size_t j = 0; j < max_work_item_dimensions; ++j) {
                        max_work_item_sizes << " " << (device_m[i].getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>())[j];
                    }
                    std::cerr << "  device[" << i << "]:\n"
                              << "    type: " << device_type << "\n"
                              << "    vendor: " << device_m[i].getInfo<CL_DEVICE_VENDOR>() << "\n"
                              << "    name: " << device_m[i].getInfo<CL_DEVICE_NAME>() << "\n"
                              << "    version: " << device_m[i].getInfo<CL_DEVICE_VERSION>() << "\n"
                              << "    profile: " << device_m[i].getInfo<CL_DEVICE_PROFILE>() << "\n"
                              << "    extensions: " << device_m[i].getInfo<CL_DEVICE_EXTENSIONS>() << "\n"
                              << "    address bits: " << device_m[i].getInfo<CL_DEVICE_ADDRESS_BITS>() << "\n"
                              << "    endian little: " << device_m[i].getInfo<CL_DEVICE_ENDIAN_LITTLE>() << "\n"
                              << "    error correction support: " << device_m[i].getInfo<CL_DEVICE_ERROR_CORRECTION_SUPPORT>() << "\n"
                              << "    max compute units: " << device_m[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << "\n"
                              << "    max work group size: " << device_m[i].getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << "\n"
                              << "    max work item dimensions: " << max_work_item_dimensions << "\n"
                              << "    max work item sizes: " << max_work_item_sizes.str() << "\n"
                              << "    global mem size: " << device_m[i].getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << "\n"
                              << "    global mem cache size: " << device_m[i].getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>() << "\n"
                              << "    global mem cache type: " << device_m[i].getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>() << "\n"
                              << "    global mem cacheline size: " << device_m[i].getInfo<CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE>() << "\n"
                              << "    local mem size: " << device_m[i].getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << "\n"
                              << "    timer resolution (ns): " << device_m[i].getInfo<CL_DEVICE_PROFILING_TIMER_RESOLUTION>() << "\n"
                              << "    preferred vector width char: " << device_m[i].getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR>() << "\n"
                              << "    preferred vector width short: " << device_m[i].getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT>() << "\n"
                              << "    preferred vector width int: " << device_m[i].getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT>() << "\n"
                              << "    preferred vector width long: " << device_m[i].getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG>() << "\n"
                              << "    preferred vector width float: " << device_m[i].getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT>() << "\n"
                              << "    preferred vector width double: " << device_m[i].getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE>() << "\n"
                              << "    driver version: " << device_m[i].getInfo<CL_DRIVER_VERSION>() << "\n";

                }
            }

            cl_uint properties = 0;
            if (profile_m) {
                properties |= CL_QUEUE_PROFILING_ENABLE;
            }
            for (size_t i = 0; i < device_m.size(); ++i) {
                queue_m.push_back(cl::CommandQueue(context_m, device_m[i], properties, &rc));
            }
		}

    ~AppBase()
		{
		}

	virtual void build_program(std::string const & device_program_text = "")
		{
			if (device_program_text.size()) {
				device_program_text_m = device_program_text;
			} else {
				device_program_text_m = get_device_program_text();
			}
			
			// Build the program
			cl::Program::Sources source(1, std::make_pair(device_program_text_m.c_str(), device_program_text_m.length()));
			program_m = cl::Program(context_m, source);
			try {
				program_m.build(device_m);
			}
			catch (cl::Error error) {
				std::cerr << "ERROR: Build failed\n";
				std::cerr << "Device program text:\n" << device_program_text_m << "\n\n\n";
				for (size_t i = 0; i < device_m.size(); ++i) {
					std::string build_log = program_m.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device_m[i]);
					std::cerr << "Build log for device[" << i << "]:\n" << build_log << "\n";
				}
				throw;
			}			
		}
		
    static const char *opencl_error_string(cl_int err)
	{
		switch (err) {
			case CL_SUCCESS:                          return "Success!";
			case CL_DEVICE_NOT_FOUND:                 return "Device not found.";
			case CL_DEVICE_NOT_AVAILABLE:             return "Device not available";
			case CL_COMPILER_NOT_AVAILABLE:           return "Compiler not available";
			case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return "Memory object allocation failure";
			case CL_OUT_OF_RESOURCES:                 return "Out of resources";
			case CL_OUT_OF_HOST_MEMORY:               return "Out of host memory";
			case CL_PROFILING_INFO_NOT_AVAILABLE:     return "Profiling information not available";
			case CL_MEM_COPY_OVERLAP:                 return "Memory copy overlap";
			case CL_IMAGE_FORMAT_MISMATCH:            return "Image format mismatch";
			case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return "Image format not supported";
			case CL_BUILD_PROGRAM_FAILURE:            return "Program build failure";
			case CL_MAP_FAILURE:                      return "Map failure";
			case CL_INVALID_VALUE:                    return "Invalid value";
			case CL_INVALID_DEVICE_TYPE:              return "Invalid device type";
			case CL_INVALID_PLATFORM:                 return "Invalid platform";
			case CL_INVALID_DEVICE:                   return "Invalid device";
			case CL_INVALID_CONTEXT:                  return "Invalid context";
			case CL_INVALID_QUEUE_PROPERTIES:         return "Invalid queue properties";
			case CL_INVALID_COMMAND_QUEUE:            return "Invalid command queue";
			case CL_INVALID_HOST_PTR:                 return "Invalid host pointer";
			case CL_INVALID_MEM_OBJECT:               return "Invalid memory object";
			case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return "Invalid image format descriptor";
			case CL_INVALID_IMAGE_SIZE:               return "Invalid image size";
			case CL_INVALID_SAMPLER:                  return "Invalid sampler";
			case CL_INVALID_BINARY:                   return "Invalid binary";
			case CL_INVALID_BUILD_OPTIONS:            return "Invalid build options";
			case CL_INVALID_PROGRAM:                  return "Invalid program";
			case CL_INVALID_PROGRAM_EXECUTABLE:       return "Invalid program executable";
			case CL_INVALID_KERNEL_NAME:              return "Invalid kernel name";
			case CL_INVALID_KERNEL_DEFINITION:        return "Invalid kernel definition";
			case CL_INVALID_KERNEL:                   return "Invalid kernel";
			case CL_INVALID_ARG_INDEX:                return "Invalid argument index";
			case CL_INVALID_ARG_VALUE:                return "Invalid argument value";
			case CL_INVALID_ARG_SIZE:                 return "Invalid argument size";
			case CL_INVALID_KERNEL_ARGS:              return "Invalid kernel arguments";
			case CL_INVALID_WORK_DIMENSION:           return "Invalid work dimension";
			case CL_INVALID_WORK_GROUP_SIZE:          return "Invalid work group size";
			case CL_INVALID_WORK_ITEM_SIZE:           return "Invalid work item size";
			case CL_INVALID_GLOBAL_OFFSET:            return "Invalid global offset";
			case CL_INVALID_EVENT_WAIT_LIST:          return "Invalid event wait list";
			case CL_INVALID_EVENT:                    return "Invalid event";
			case CL_INVALID_OPERATION:                return "Invalid operation";
			case CL_INVALID_GL_OBJECT:                return "Invalid OpenGL object";
			case CL_INVALID_BUFFER_SIZE:              return "Invalid buffer size";
			case CL_INVALID_MIP_LEVEL:                return "Invalid mip-map level";
			default:                                  return "Unknown";
		}
		return "Unknown error";
    }

	// find the most capable device based on type, and then number of compute units
	int get_most_capable_device()
	{
		int device_id = -1;
		std::vector<unsigned int> device_type;
		device_type.push_back(CL_DEVICE_TYPE_GPU);
		device_type.push_back(CL_DEVICE_TYPE_ACCELERATOR);
		device_type.push_back(CL_DEVICE_TYPE_CPU);
		
		std::vector<unsigned int>::iterator type;
		for (type = device_type.begin(); type != device_type.end(); ++type) {
			unsigned int max_compute_units = 0;
			for (size_t i = 0; i < device_m.size(); ++i) {
				if (*type == device_m[i].getInfo<CL_DEVICE_TYPE>()) {
					if (max_compute_units < device_m[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>()) {
						max_compute_units = device_m[i].getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
						device_id = i;
					}
				}
			}
			if (device_id >= 0) {
				break;
			}
		}
		
		return device_id;
	}
	
    virtual void host_run() = 0;

protected:
	
	virtual std::string const & get_device_program_text() = 0;

    cl::Context context_m;
    cl::Program program_m;
    int debug_m;
    int profile_m;
    int verbose_m;
    std::string device_program_text_m;
    std::vector<cl::CommandQueue> queue_m;
    std::vector<cl::Device> device_m;
    std::vector<cl::Platform> platform_m;
	
};

#endif
