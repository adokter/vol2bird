#include <torch/script.h> 
#include <iostream>
#include <memory>

using namespace std;

extern "C" int run_mistnet(float* tensor_in, float** tensor_out, const char* model_path, int tensor_size) {
    
        // ***************************************************************************
        // *************************                           ***********************
        // ************************* the code to use the model ***********************
        // *************************                           ***********************
        // ***************************************************************************
        
        torch::jit::script::Module module;
        try {
            module = torch::jit::load(model_path);
        }
        catch (const c10::Error& e) {
            std::cerr << "\nError: failed to load MistNet model from file " << model_path << "\n";
            return -1;
        }

        // if you already have a 1d floating point array that is the tensor of size 15 x 608 x 608
        // pointed by a pointer (float*) tensor_in, you can convert it to a torch tensor by:
        at::Tensor inputs = torch::from_blob(tensor_in, {1, 15, 608, 608}, at::kFloat);

        std::vector<torch::jit::IValue> inputs_;
        inputs_.push_back(inputs);

        at::Tensor output = module.forward(inputs_).toTensor();

        float *output_array = output.data<float>();
        
        // copy result
        for (int i=0; i<3*5*608*608; i++) (*tensor_out)[i] = output_array[i];
        
        return 0;
}

