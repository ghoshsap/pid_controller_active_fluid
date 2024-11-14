#include <fstream>
#include <sstream> 
#include <vector>
#include <iostream>
#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <string>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <cuda_runtime_api.h>
#include <driver_types.h>
#include <ostream>
#include <chrono>
#include <ctime>
#include "cuPSS/inc/cupss.h"

#ifdef WITHCUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#define NX 128
#define NY 128 
#define NSTEPS 20000000

int main(int argc, char **argv)
{


    float kp = std::atof(argv[1]);   // Proportional gain
    float ki = std::atof(argv[2]);  // Integral gain
    float dt = 0.0001f; // Time step 
    float tc = std::atof(argv[3]);    // Time constant (tc)
    float err_int = 0.0f;        // Integral of error
    float v_set = 2.0f;          // Desired velocity
    // float u_mean = 0.0f;         // Average velocity (to be computed later)
    int steps_per_l = 1;       // Interval for updating the controller
    int count = 0;

    float alpha;
    float alpha_prev;

    std::vector<float> u_mean_list;
    std::vector<float> alpha_mean_list;

    evolver system(1, NX, NY, 1.0f, 1.0f, dt, 10000);

    system.createField("Qxx", true);
    system.createField("Qxy", true);
    system.createField("alpha", false);
    system.createField("iqxQxx", false);
    system.createField("iqxQxy", false);
    system.createField("iqyQxx", false);
    system.createField("iqyQxy", false);
    system.createField("sigxx", false);
    system.createField("sigxy", false);
    system.createField("vx", false);
    system.createField("vy", false);
    system.createField("w", false);
    system.createField("Q2", false);

    system.addParameter("a2", -1.0f);
    system.addParameter("a4", 1.0f);
    system.addParameter("kQ", 4.0f);
    system.addParameter("eta", 1.0f);
    system.addParameter("gamma", 1.0f);
    system.addParameter("lambda", 1.0f);
    system.addParameter("fric", 0.01);



    system.addEquation("dt Qxx + (a2*1/gamma + kQ*1/gamma*q^2)*Qxx = -a4*1/gamma*Q2*Qxx - vx*iqxQxx - vy*iqyQxx + lambda*iqx*vx - 2*Qxy*w");
    system.addEquation("dt Qxy + (a2*1/gamma + kQ*1/gamma*q^2)*Qxy = -a4*1/gamma*Q2*Qxy - vx*iqxQxy - vy*iqyQxy + 0.5*lambda*iqx*vy + 0.5*lambda*iqy*vx + 2*Qxx*w");
    system.addEquation("iqxQxx = iqx*Qxx");
    system.addEquation("iqxQxy = iqx*Qxy");
    system.addEquation("iqyQxx = iqy*Qxx");
    system.addEquation("iqyQxy = iqy*Qxy");

    system.addEquation("sigxx = alpha * Qxx");
    system.addEquation("sigxy = alpha * Qxy");
    system.addEquation("w = 0.5*iqx*vy-0.5*iqy*vx");
    system.addEquation("Q2 = Qxx^2 + Qxy^2");

    system.addEquation("vx * (fric + eta*q^2) = (iqx + iqx^3*1/q^2 - iqx*iqy^2*1/q^2) * sigxx + (iqy + iqx^2* iqy*1/q^2 + iqx^2*iqy*1/q^2) * sigxy");
    system.addEquation("vy * (fric + eta*q^2) = (iqx + iqx*iqy^2*1/q^2 + iqx*iqy^2*1/q^2) * sigxy + (-iqy - iqy^3*1/q^2 + iqx^2*iqy*1/q^2) * sigxx");

    system.printInformation();

    // system.addNoise("Qxx", "1.5*q^2");
    // system.addNoise("Qxy", "1.5*q^2");    

    // Random initial state
    std::srand(1324);
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            int index = j * NX + i;
            system.fields[0]->real_array[index].x = -0.0f + 0.001f * (float)(rand()%200-100);
            system.fields[0]->real_array[index].y = 0.0f;
        }
    }

    system.prepareProblem();

    system.fields[0]->outputToFile = true;
    for (int i = 0; i < system.fields.size(); i++)
    {
        system.fields[i]->outputToFile = false;
    }
    system.setOutputField("Q2", true);
    system.setOutputField("Qxx", true);
    system.setOutputField("Qxy", true);
    system.setOutputField("vx", true);
    system.setOutputField("vy", true);
    system.setOutputField("alpha", true);


    float2 alpha_temp[NX*NY];

    float2 vx_temp[NX*NY];
    float2 vy_temp[NX*NY];

    for (int i = 0; i < NX*NY; i++)
    {        
        alpha_temp[i].x = 0.0f;
        alpha_temp[i].y = 0.0f;
    }


    int steps = NSTEPS;
    int check = steps/100;
    if (check < 1) check = 1;


    auto start = std::chrono::system_clock::now();
    int frame_count = 0;

    for (int i = 0; i < steps; i++) {

        // float error = v_set - u_mean;     
        // err_int += error;   

        if (i % steps_per_l == 0) {

            cudaMemcpy(vx_temp, system.fields[9]->real_array_d, NX*NY*sizeof(float2), cudaMemcpyDeviceToHost);
            cudaMemcpy(vy_temp, system.fields[10]->real_array_d, NX*NY*sizeof(float2), cudaMemcpyDeviceToHost);


            //////////////////////$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$/////////////////////////////

            // Measure u_mean (average velocity)
            double u_mean = 0.0f;
            for (int j = 0; j < NY; j++) {
                for (int k = 0; k < NX; k++) {
                    int index = j * NX + k;
                    float vx = vx_temp[index].x;  
                    float vy = vy_temp[index].x;  
                    u_mean += sqrt(vx * vx + vy * vy);  
                }
            }
            u_mean /= (NX * NY);  

            //////////////////////$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$/////////////////////////////

            alpha_mean_list.push_back(alpha);
            u_mean_list.push_back(u_mean); // half domain

       
            // ###### PI CONTROLLER #######
            if (count == 0) {
                
                alpha = -0.2f;
                alpha_prev = alpha;
            }

            if (count % steps_per_l == 0 && count > 0) {

                float error = v_set - u_mean;     
                err_int += error;    

                float pi_control = -kp * error - ki * dt * steps_per_l* err_int ;

                // alpha = pi_control;
                alpha = alpha_prev + (dt*steps_per_l ) * (-alpha_prev/(tc) + pi_control/tc);

                float alpha_cut = std::min(alpha, 0.0f);
                alpha = std::max(alpha, alpha_cut);

                alpha_prev = alpha;
            }


            int index = 0;
            for (int i = 0; i < NX; i++)
            {
                for (int j = 0; j < NY; j++)
                {
                    index = j*NX + i;
                    alpha_temp[index].x = alpha;

                }
            }

            // cudaMemcpy(destination pointer, origin pointer, size, direction);
            cudaMemcpy(system.fields[2]->real_array_d, alpha_temp, NX*NY*sizeof(float2), cudaMemcpyHostToDevice);
            // cudaMemcpy(system.fields[2]->real_array, alpha_temp, NX*NY*sizeof(float2), cudaMemcpyHostToHost);
            system.fields[2]->toComp();

        }

        system.advanceTime();

        count++;

        if (i % (steps / 100) == 0) {
            std::cout << "Progress: " << i / (steps / 100) << "%\r";
            std::cout.flush();
        }
    }

    std::cout << "Simulation finished." << std::endl;
 
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::ostringstream ss;
    ss.precision(2);  // Set precision to one decimal place
    ss << std::fixed << "output_means_kp_" << kp << "_ki_" << ki << "_tc_" << tc << ".txt";

    std::string filename = ss.str();  // Convert the stream to a string
    
    std::ofstream outfile(filename);
    for (size_t i = 0; i < u_mean_list.size(); i++) {
        outfile << "u_mean = " << u_mean_list[i] << ", alpha_mean = " << alpha_mean_list[i] << "\n";
    }
    outfile.close();
 
 
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s"
              << std::endl;
    return 0;
}
