#======================================================================
    Development script for mgm.jl
    
    This script trains an autoencoder neural network to find a good
    parametrization of the pressure distribution cp(x)

    Author: Andrea Pavan
    Project: AirfoilInverseDesign.jl package
    License: MIT
    Date: 30/12/2023
======================================================================#
using Flux;
#using JLD2;
using Plots;

println("AirfoilInverseDesign - Training neural network");
datasetfile = joinpath(@__DIR__,"02_dataset_x20.csv");
#batchsize = 65;

#count the number of samples in the CSV dataset
N = 0;                              #total number of samples in the dataset
for line in eachline(datasetfile)
    if line[1]=='"'
        global N += 1;
    end
end
println("Found ",N, " samples in the dataset");
Ntrain = floor(Int,0.8*N);          #number of training samples
Nvalidation = N-Ntrain;             #number of validation samples
#drop some samples in order to have a constant batch size
#Ntrain -= Ntrain%batchsize;
#Nvalidation -= Nvalidation%batchsize;


#import data
ytrain = zeros(Float32,(20,Ntrain));
yvalidation = zeros(Float32,(20,Nvalidation));
i = 1;
#read dataset from CSV
for line in eachline(datasetfile)
    if line[1]=='"'
        currententry = parse.(Float32,split(replace(line,"["=>"","]"=>""," "=>""), ",")[2:end]);
        if i<=Ntrain
            ytrain[:,i] = currententry;
        elseif Ntrain<i<=Ntrain+Nvalidation
            yvalidation[:,i-Ntrain] = currententry;
        end
        global i += 1;
    end
end
println("Dataset parsed successfully");


#define neural network
println("Defining a MLP neural network:");
myneuralnet = Chain(
    Dense(20 => 12),
    Dense(12 => 8, tanh),
    Dense(8 => 12, tanh),
    Dense(12 => 20)
);
#OR resume existing neural network from file
#myneuralnet = JLD2.load("dev/02_myneuralnet_training_checkpoint.jld2","myneuralnet");
display(myneuralnet);
lossfun(y,yexact) = Flux.mse(y,yexact);
#lossfun(y,yexact) = Flux.mae(y,yexact);


#initial evaluation
println("Initial inference:");
ynn = myneuralnet(ytrain);
println("  loss = ",lossfun(ynn,ytrain));


#training
println("Training:");
epochs = 20000;                                  #number of epochs
trainingloss = zeros(Float64,epochs);           #training loss history over epochs
validationloss = zeros(Float64,epochs);         #validation loss history over epochs
#trainingdata = Flux.DataLoader((xtrain,ytrain),batchsize=batchsize);        #setting up training data for batch use
optimizer = Flux.setup(Adam(), myneuralnet);                                #Adam optimizer without regularization
@time for epoch in 1:epochs
    #single batch training
    (currentloss,currentgradient) = Flux.withgradient(m -> lossfun(m(ytrain),ytrain), myneuralnet);
    Flux.update!(optimizer, myneuralnet, currentgradient[1]);

    #minibatch training
    #=currentloss = 0;
    for (xbatch,ybatch) in trainingdata
        (currentloss,currentgradient) = Flux.withgradient(m -> lossfun(m(ybatch),ybatch), myneuralnet);
        Flux.update!(optimizer, myneuralnet, currentgradient[1]);
    end=#

    #evaluate performance of the current model
    trainingloss[epoch] = currentloss;
    ynnval = myneuralnet(yvalidation);
    validationloss[epoch] = lossfun(ynnval,yvalidation);
    if epoch%100 == 0
        println("  epoch ",epoch," - training loss = ",rpad(round(trainingloss[epoch],digits=6),8,'0')," - validation loss = ",rpad(round(validationloss[epoch],digits=6),8,'0'));
    end

    #save current model if it outperforms the best model so far
    #=if epoch>5000 && validationloss[epoch]<=minimum(validationloss[1:epoch-1])
        jldsave("dev/02_myneuralnet_training_checkpoint.jld2"; myneuralnet);
    end=#
end
println("Training completed");

#training convergence plot
plt2 = plot(1:epochs, trainingloss, color=:blue, linewidth=3, label="Training loss",
    title="Training & validation loss",
    xlabel="Epoch",
    ylabel="Loss function"
);
plot!(plt2, 1:epochs, validationloss, color=:orange, linewidth=3, label="Validation loss");
display(plt2);


#save model
jldsave("dev/02_myneuralnet.jld2"; myneuralnet);
#mv("dev/02_myneuralnet_training_checkpoint.jld2","dev/02_myneuralnet.jld2",force=true);
#println("Model saved to file");


#display neural network parameters
println("myneuralnet parameters:");
for layer in myneuralnet.layers
    display(layer);
    println("> W = ",string(string(layer.weight)),"\n");
    println("> b = ",string(layer.bias),"\n\n");
end
