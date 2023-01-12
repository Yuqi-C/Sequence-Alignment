# Sequence-Alignment

## Files
1. basic.cc implements the basic dynamic programming solution to the problem.
2. efficient.cc provides the memory efficient solution.
3. 2 Shell files ‘basic.sh’ and ‘ efficient.sh’ are the commands to compile and run two cpp files.
4. Folder 'datapoints' includes .txt files used to generate string based on 

## Input string Generator
The input to the program would be a text file containing the following information:
1. First base string (𝑠1)
2. Next 𝑗 lines consist of indices after which the copy of the previous string needs to be inserted in the cumulative string. (eg given below)
3. Second base string (𝑠2)
4. Next 𝑘 lines consist of indices after which the copy of the previous

## Values for Delta and Alphas

Values for α’s are as follows. δ𝑒 is equal to 30.

|      |   A   |   C   |   G   |   T   |
|:----:|:-----:|:-----:|:-----:|:-----:|
|  **A**   |   0   |   110   |   48   |   94   |
|  **C**   |   110  |   0   |   118   |   48   |
|  **G**   |   48   |   118   |   0   |   110   |
|  **T**   |  94   |   48   |   110   |   0   |
