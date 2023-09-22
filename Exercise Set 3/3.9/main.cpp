#include <chrono>
#include "hoch_helper.h"
int main () {
    auto start_time = chrono::high_resolution_clock::now();
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end_time - start_time;
    
    FILE* my_file = fopen("test.txt","w");
    
    // test quat_mult
    double q1[4] = {0.5,0.5,0.5,0.5};
    double q2[4];
    double q3[4];
    double euler[3] = {0.0,0.0,0.0};
    array_copy(q1,q3,4);
    
    start_time = chrono::high_resolution_clock::now();
    for (int i=1; i<11; i++){
        for (int j=1; j<100001; j++){
            array_copy(q3,q2,4);
            for (int k=0; k<100; k++){
                quat_mult(q1,q2,q3);
            }
        }
    fprintf(my_file,"%5d%25.11E%25.11E%25.11E%25.11E\n",i,q3[0],q3[1],q3[2],q3[3]);
    }
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    printf("quat_mult time total [sec]: %f\n",duration.count());

    // test euler_to_quat
    start_time = chrono::high_resolution_clock::now();
    for(int i=-180; i< 182; i+=2){
        euler[0] = float(i)*pi/180.0;
        for (int j=-90; j<92; j+=2){
            euler[1] = float(j)*pi/180.0;
            for (int k=0; k<362; k+=2){
                euler[2] = float(k)*pi/180.0;
                for (int m=0; m<10; m++){
                    euler_to_quat(euler,q1);
                }
            }
        }
        fprintf(my_file,"%5d%25.11E%25.11E%25.11E%25.11E\n",i,q1[0],q1[1],q1[2],q1[3]);
    }
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    printf("euler_to_quat time total [sec]: %f\n",duration.count());

    // test quat_to_euler
    start_time = chrono::high_resolution_clock::now();
    for (int i=4; i<102; i+=4){
        for (int j=4; j<102; j+=4){
            for (int k=4; k<102; k+=2){
                for (int m=4; m<102; m+=2){
                    q1[0] = 0.01*float(i);
                    q1[1] = 0.01*float(j);
                    q1[2] = 0.01*float(k);
                    q1[3] = 0.01*float(m);
                    quat_norm(q1);
                    for(int n=0; n<10; n++){
                        quat_to_euler(q1,euler);
                    }
                }
            }
        }
        fprintf(my_file,"%5d%25.11E%25.11E%25.11E\n",i,euler[0]*180.0/pi, euler[1]*180.0/pi, euler[2]*180.0/pi);
    }
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    printf("quat_to_euler time total [sec]: %f\n",duration.count());

    // test body_to_fixed
    double vec1[3] = {0.,0.,1.};
    double vec2[3];
    start_time = chrono::high_resolution_clock::now();
    for (int i=1; i<102; i+=4){
        for (int j=1; j<102; j++){
            for (int k=1; k<102; k++){
                for (int m=1; m<102; m++){
                    q1[0] = 0.01*float(i);
                    q1[1] = 0.01*float(j);
                    q1[2] = 0.01*float(k);
                    q1[3] = 0.01*float(m);
                    body_to_fixed(vec1,q1,vec2);
                }
            }
        }
        fprintf(my_file,"%5d%25.11E%25.11E%25.11E\n",i,vec2[0],vec2[1],vec2[2]);
    }
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    printf("body_to_fixed time total [sec]: %f\n",duration.count());

    // test fixed_to_body
    vec1[0] = 0.0;
    vec1[1] = 0.0;
    vec1[2] = 1.0;
    start_time = chrono::high_resolution_clock::now();
    for (int i=1; i<102; i+=4){
        for (int j=1; j<102; j++){
            for (int k=1; k<102; k++){
                for (int m=1; m<102; m++){
                    q1[0] = 0.01*float(i);
                    q1[1] = 0.01*float(j);
                    q1[2] = 0.01*float(k);
                    q1[3] = 0.01*float(m);
                    fixed_to_body(vec1,q1,vec2);
                }
            }
        }
        fprintf(my_file,"%5d%25.11E%25.11E%25.11E\n",i,vec2[0],vec2[1],vec2[2]);
    }
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    printf("fixed_to_body time total [sec]: %f\n",duration.count());

    fclose(my_file);
    return 0;
}