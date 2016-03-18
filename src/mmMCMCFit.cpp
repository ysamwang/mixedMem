#include "mmMCMCFit.h"

void mixedMemMcmcC(mm_model_mcmc model, int burnIn, int samples, int thin, int print, Rcpp::List fileNames, int newFiles,
                   double omega, double eta, NumericVector which_write)
{
    if(newFiles) {
        clearFiles(fileNames, which_write);
    }

    int t,s;

    for(t = 0; t < burnIn; t++) {
        if(t% print == 0) {
            Rcout<<"Burn-in: "<<t<<std::endl;
        }

        //Data Augmentation
        if(model.getExtended() == 1) {
            imputeStayers(model);
        }

        imputeLatent(model);


        // posterior draws of model parameters
        sampleTheta(model) ;

        sampleKsi(model, eta);
        // Draw dirichlet parameter
        sampleAlpha(model, omega);


        // Draw Extended
        if(model.getExtended() == 1) {
            sampleP(model);
        }
    }

    Rcout<<"=== Done with Burn In ==="<<std::endl;
    for(t = 0; t < samples; t++) {
        for(s = 0; s < thin; s++) {

            //Data Augmentation
            if(model.getExtended() == 1) {
                imputeStayers(model);
            }
            imputeLatent(model);

            // posterior draws of model parameters
            sampleTheta(model);

            sampleKsi(model, eta);

            // Draw dirichlet parameter
            sampleAlpha(model, omega);

            // Sample P
            if(model.getExtended() == 1) {
                sampleP(model);
            }

        }

        if(t% print == 0) {
            Rcout<<"Sample : "<<t<<std::endl;
        }

        writeState(model, fileNames, which_write);
    }

    Rcout<<"Sampling Completed!"<<std::endl;
}


void writeState(mm_model_mcmc model, Rcpp::List fileNames, NumericVector which_write)
{
    int i, j, k, r, v;
    std::ofstream myfile;
    std::string file;


    //Write Theta
    if(which_write[0] == 1) {
        file = as<std::string> (fileNames[0]);
        myfile.open (file.c_str(),std::ios::app);
        for(v = 0; v < model.getMaxV(); v++) {
            for(k = 0; k < model.getK(); k++) {
                for(j = 0; j < model.getJ(); j++) {
                    if((j + k + v) > 0) {
                        myfile<< ", ";
                    }
                    if(v < model.getV(j)) {
                        myfile << model.getTheta(j, k, v);
                    } else {
                        myfile << 0;
                    }
                }
            }
        }
        myfile<<std::endl;
        myfile.close();
    }

// Write Alpha
    if(which_write[1] == 1) {
        file = as<std::string> (fileNames[1]);
        myfile.open (file.c_str(),std::ios::app);
        myfile << model.getAlpha();
        myfile << std::endl;
        myfile.close();
    }

    //Write ksi
    if(which_write[2] == 1) {
        //Write ksi
        file = as<std::string> (fileNames[2]);
        myfile.open(file.c_str(),std::ios::app);
        for(k = 0; k <model.getK(); k++) {
            if(k > 0) {
                myfile<< ", ";
            }
            myfile<< model.getKsi(k);
        }
        myfile<<std::endl;
        myfile.close();

    }


    //Write Lambda
    if(which_write[3] == 1) {
        file = as<std::string> (fileNames[3]);
        myfile.open (file.c_str() ,std::ios::app);
        for(k = 0; k < model.getK(); k++) {
            for(i = 0; i < model.getT(); i++) {
                if((i+k) > 0) {
                    myfile<< ", ";
                }
                myfile << model.getLambda(i,k);
            }
        }
        myfile << std::endl;
        myfile.close();
    }

    //write Z
    if(which_write[4] == 1) {

        file = as<std::string> (fileNames[4]);
        myfile.open (file.c_str(),std::ios::app);
        for(r = 0; r < model.getMaxR(); r++) {
            for(j = 0; j < model.getJ(); j++) {
                for(i = 0; i < model.getT(); i++) {
                    if((i + j + r) > 0) {
                        myfile<< ", ";
                    }
                    if(r < model.getR(j)) {
                        myfile << model.getZ(i,j,r);
                    } else {
                        myfile << -1 ;
                    }
                }
            }
        }
        myfile<<std::endl;
        myfile.close();

    }


    // Write P
    if(which_write[5] == 1) {
        file = as<std::string> (fileNames[5]);
        myfile.open (file.c_str(),std::ios::app);
        for(k = 0; k < (model.getS() + 1); k++) {
            if(k > 0) {
                myfile<< ", ";
            }
            myfile << model.getP(k);
        }
        myfile<<std::endl;
        myfile.close();
    }

    // Write stayerStatus
    if(which_write[6] == 1) {
        file = as<std::string> (fileNames[6]);
        myfile.open (file.c_str(),std::ios::app);
        for(i = 0; i < model.getT() ; i ++) {
            if(i > 0) {
                myfile << ", ";
            }
            myfile << model.getStayerStatus(i);
        }
        myfile<<std::endl;
        myfile.close();

    }
}

void clearFiles(Rcpp::List fileNames, NumericVector which_write)
{
    std::ofstream myfile;
    std::string file;


    if(which_write[0] == 1) {
        file = as<std::string> (fileNames[0]);
        myfile.open(file.c_str());
        myfile.close();
    }

    if(which_write[1] == 1) {
        file = as<std::string> (fileNames[1]);
        myfile.open(file.c_str());
        myfile.close();
    }

    if(which_write[2] == 1) {
        file = as<std::string> (fileNames[2]);
        myfile.open(file.c_str());
        myfile.close();
    }

    if(which_write[3] == 1) {
        file = as<std::string> (fileNames[3]);
        myfile.open(file.c_str());
        myfile.close();
    }

    if(which_write[4] == 1) {
        file = as<std::string> (fileNames[4]);
        myfile.open(file.c_str());
        myfile.close();
    }

    if(which_write[5] == 1) {
        file = as<std::string> (fileNames[5]);
        myfile.open(file.c_str());
        myfile.close();
    }

    if(which_write[6] == 1) {
        file = as<std::string> (fileNames[6]);
        myfile.open(file.c_str());
        myfile.close();
    }
}
