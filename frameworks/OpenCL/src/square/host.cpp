// Host code...

#include "square/app.hpp"

void App::host_run()
{
    int rc;

    // host data for tests
    size_t const size = 1 << 25;
    std::vector<float> idata;
    std::vector<float> odata;

    for (size_t i = 0; i < size; ++i) {
        idata.push_back((i));
    }
    odata.resize(idata.size());

    // Events for timing
    cl::Event event1, event2, event3;

    // Allocate device memory
    cl::Buffer ibuf(context_m, CL_MEM_READ_ONLY, sizeof(float) * idata.size());
    cl::Buffer obuf(context_m, CL_MEM_WRITE_ONLY, sizeof(float) * odata.size());

    // Create kernel
    cl::Kernel kernel(program_m, "square", &rc);
    rc = kernel.setArg(0, ibuf);
    rc = kernel.setArg(1, obuf);
    rc = kernel.setArg(2, size);

    // Select a device
    int device_id;
    if (device_list_m.size()) {
        // use first device in the device list
        device_id = device_list_m[0];
    } else {
        device_id = get_most_capable_device();
    }
    if (verbose_m) {
        std::cerr << "[Device selection]\n  using device[" << device_id << "]:\n";
    }

    // Copy input
    queue_m[device_id].enqueueWriteBuffer(ibuf,
                                          CL_TRUE,
                                          0,
                                          idata.size() * sizeof(float),
                                          &idata[0],
                                          NULL,
                                          &event1);

    // Run kernel
    size_t work_group_size;
    work_group_size = kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device_m[device_id]);
    if (verbose_m) {
        std::cerr << "  work group size = " << work_group_size << "\n";
    }
    queue_m[device_id].enqueueNDRangeKernel(kernel,
                                            cl::NullRange,
                                            cl::NDRange(size),
                                            cl::NDRange(work_group_size),
                                            NULL,
                                            &event2);
    // Copy output
    queue_m[device_id].enqueueReadBuffer(obuf,
                                         CL_TRUE,
                                         0,
                                         odata.size() * sizeof(float),
                                         &odata[0],
                                         NULL,
                                         &event3);

    // Wait for completion
    queue_m[device_id].finish();

    // Timings
    if (profile_m) {
        cl_ulong start, end;
        float t; // execution time in milliseconds

        std::cerr << "[Timing]\n";
        start = event1.getProfilingInfo<CL_PROFILING_COMMAND_START>();
        end = event1.getProfilingInfo<CL_PROFILING_COMMAND_END>();
        t = (end - start) * 1.0e-6f;
        std::cerr << "  write execution time = " << t << " ms\n";
        start = event2.getProfilingInfo<CL_PROFILING_COMMAND_START>();
        end = event2.getProfilingInfo<CL_PROFILING_COMMAND_END>();
        t = (end - start) * 1.0e-6f;
        std::cerr << "  kernel execution time = " << t << " ms\n";
        start = event3.getProfilingInfo<CL_PROFILING_COMMAND_START>();
        end = event3.getProfilingInfo<CL_PROFILING_COMMAND_END>();
        t = (end - start) * 1.0e-6f;
        std::cerr << "  read execution time = " << t << " ms\n";
    }

    // Check results
    std::cout << "[Results]\n";
    int nerror = 0;
    for (size_t i = 0; i < idata.size(); ++i) {
        float result = idata[i] * idata[i];
        if (result != odata[i]) {
            nerror++;
            std::cout << "  error: idata[" << i << "]^2 ("
                      << result
                      << ") != odata["
                      << i
                      << "] ("
                      << odata[i]
                      << ")\n";
        }
    }
    std::cerr << "  found " << nerror << " error(s) out of " << idata.size() << "\n";
}
