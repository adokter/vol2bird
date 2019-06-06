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
                                array[i][j][k] = data +  i * dim4 * dim3 * dim2 + j * dim4 * dim3 + k * dim4; 
                        }
                }
        }

        return array;
}

int run_mistnet(float* tensor_in, float* tensor_out, const char* model_path) {
        if (2 != 2) {
                std::cerr << "usage: example-app <path-to-exported-script-module>\n";
                return -1;
        }

	// load the scan from the hard drives
        
        /*
        cnpy::NpyArray dz = cnpy::npz_load("../dz/KBGM20030901_034952.npz", "rdata");
        cnpy::NpyArray vr = cnpy::npz_load("../vr/KBGM20030901_034952.npz", "rdata");
        cnpy::NpyArray sw = cnpy::npz_load("../sw/KBGM20030901_034952.npz", "rdata");

        cnpy::NpyArray dz = cnpy::npz_load("../dz/KMOB20080901_040557_V03.npz", "rdata");
        cnpy::NpyArray vr = cnpy::npz_load("../vr/KMOB20080901_040557_V03.npz", "rdata");
        cnpy::NpyArray sw = cnpy::npz_load("../sw/KMOB20080901_040557_V03.npz", "rdata");


        float* dz_ = dz.data<float>();
        float* vr_ = vr.data<float>();
        float* sw_ = sw.data<float>();

	
        // convert the 1d array of floating points into a torch tensor
        at::Tensor tensor_dz = torch::from_blob(dz_, {1, 5, 608, 608}, at::kFloat);
        at::Tensor tensor_vr = torch::from_blob(vr_, {1, 5, 608, 608}, at::kFloat);
        at::Tensor tensor_sw = torch::from_blob(sw_, {1, 5, 608, 608}, at::kFloat);

        */

        // ***************************************************************************
        // *************************                           ***********************
        // *************************                           ***********************
        // ************************* the code to use the model ***********************
        // *************************                           ***********************
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
        tensor_out=output.data<float>();

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
        float* copy_output = (float *)malloc(3 * 5 * 608 * 608 *sizeof(float));
        float* temp = copy_output;
        for (int i=0 ; i < 3 ; i++)
                for (int j=0 ; j < 5 ; j++)
                        for (int k=0 ; k < 608 ; k++)
                                for (int l=0 ; l < 608 ; l++){
                                        *temp = output_array[i][j][k][l];
                                        temp++;
                                }
        //cnpy::npy_save("./copy_output.npy", copy_output, {1, 3, 5, 608, 608}, "w");
        return 0;
}
