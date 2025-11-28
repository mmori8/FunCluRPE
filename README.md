# FunCluRPE[README.md](https://github.com/user-attachments/files/23827514/README.md)
# FunCluRPE

Clustering functional data, where each observation is represented by a curve defined over a continuous domain, remains challenging due to high dimensionality and the need for stable, data-adaptive partitioning. We propose a clustering framework based on Random Projections, which simultaneously performs dimensionality reduction and generates multiple stochastic representations of the original functions.

## Included Functions

### **FunCluRPE**
The main function. Requires either a matrix of coefficients or an `fd`-class type object (from the **fda** package) and either `d` or `G` to set the dimension for the projected space. Additional optional parameters can be set.

### **data_generation**
A function to generate example curves as in Scenario 1 of the paper.

### **run_example**
An illustrative example file demonstrating usage.

