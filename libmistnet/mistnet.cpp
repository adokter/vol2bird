#include <torch/script.h> 
#include <iostream>
#include <memory>
//#include"/usr/local/include/cnpy.h"
//#include"cnpy.h"


using namespace std;
// TODO: need to free the memory for pointers
float**** create4DTensor( float* data, int dim1, int dim2, int dim3, int dim4) {
        float ****array = (float ****)malloc(dim1 * sizeof(float***));
        for(int i=0 ; i < dim1 ; i++) {
                array[i] = (float ***) malloc(dim2 * sizeof(float**));
                for(int j=0 ; j < dim2 ; j++) {
                        array[i][j] = (float **)malloc(dim3 * sizeof(float*));
                        for (int k=0 ; k < dim3 ; k++){
                        	array[i][j][k] = (float *)malloc(dim4 * sizeof(float));
				for (int l=0 ; l < dim4 ; l++){
                                	array[i][j][k][l] = data[i * dim4 * dim3 * dim2 + j * dim4 * dim3 + k * dim4 + l]; 
				}
                        }
                }
        }

        return array;
}

extern "C" int run_mistnet(float* tensor_in, float** tensor_out, const char* model_path) {
        if (2 != 2) {
                std::cerr << "usage: example-app <path-to-exported-script-module>\n";
                return -1;
        }

        // ***************************************************************************
        // *************************                           ***********************
        // ************************* the code to use the model ***********************
        // *************************                           ***********************
        // ***************************************************************************
        std::shared_ptr<torch::jit::script::Module> module = torch::jit::load(model_path);

        assert(module != nullptr);
        
        // if you already have a 1d floating point array that is the tensor of size 15 x 608 x 608
        // pointed by a pointer (float*) tensor_in, you can convert it to a torch tensor by:
        at::Tensor inputs = torch::from_blob(tensor_in, {1, 15, 608, 608}, at::kFloat);

	//at::Tensor inputs = at::cat({tensor_dz, tensor_vr, tensor_sw}, 1);

        std::vector<torch::jit::IValue> inputs_;
        inputs_.push_back(inputs);

        at::Tensor output = module->forward(inputs_).toTensor();

        float**** output_array = create4DTensor(output.data<float>(), 3, 5, 608, 608);

        // ***************************************************************************
        // ***************************************************************************
        // ***************************************************************************
        // ***************************************************************************
        // ***************************************************************************

	//cnpy::npy_save("./output.npy", output.data<float>(), {1, 3, 5, 608, 608}, "w");
        /* sanity check1
        cout << output_array[1][0][200][200] << endl;
        cout << output_array[2][1][300][300] << endl;
        */

        // sanity check
        //float* copy_output = (float *)malloc(3 * 5 * 608 * 608 *sizeof(float));
        //float* temp = copy_output;
        for (int i=0 ; i < 3 ; i++){
                for (int j=0 ; j < 5 ; j++){
                        for (int k=0 ; k < 608 ; k++){
                                for (int l=0 ; l < 608 ; l++){
                                        (*tensor_out)[i*5*608*608 + j*608*608 +k*608 +l] = output_array[i][j][k][l];
                                }
			}
		}
	}
        //cnpy::npy_save("./copy_output.npy", copy_output, {1, 3, 5, 608, 608}, "w");
        return 0;
}
