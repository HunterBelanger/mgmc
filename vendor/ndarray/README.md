# NDArray: Numpy style arrays for C++
A common scenario which arrises in scientific and high performance computing
is the need to conduct calculations or simulations in faster, lower level
languages such as C++, and then load results into a higher level language such
as python to conduct analysis on the results. Depending on the type and amount
of data, it can be difficult and time consuming to write this data to a file in
a specific format, making sure it is consisten, then writing a second python
implementation to reliably read such data. If storing as raw text however, this
can lead to long load times, and large data files.

Currently, there are few small and complete C++ libraries to complete this task.
The xtensor library is a very good option, and I do recommend it. However, it
is a large library which is likely much more than many need, especially in 
smaller projects. If you are working on a larger project and need more
NumPy magic than this library provides, I do recommend that you check it out
[here](https://github.com/xtensor-stack/xtensor). It likely has everything
you are looking for.

I have written the NDArray template library, a minimal implementation of Numpy
style arrays for C++, as a smaller alternative to xtensor. These arrays can be
multi-dimensional, and reshaped like Numpy arrays. Indexing can be done either
with a vector, or as variadic parammeters, both using the () operator. Access
to the data using the linear index is also permitted via the [] operator.

It is also possible to load/save data from/to a ```.npy``` binary file. This
allows for fast and easy access to the data in python (as well as many other
languages). While the template container can be used to store any array of
any type, the load and save methods are only valid for the following templates:

* ```NDArray<char>```
* ```NDArray<unsigned char>```
* ```NDArray<uint16_t>```
* ```NDArray<uint32_t>```
* ```NDArray<uint64_t>```
* ```NDArray<int16_t>```
* ```NDArray<int32_t>```
* ```NDArray<int64_t>```
* ```NDArray<float>```
* ```NDArray<double>```
* ```NDArray<std::complex<float>>```
* ```NDArray<std::complex<double>>```

This is due to the fact that only certain numberic types are allowed by Numpy.
Attempting to save a NDArray templated for a different type will result in
and exception being thrown. Python and Numpy allow for the storing of raw
Python objects in ```.npy``` files, but the loading of such files into a C++
program with this library will also result in an exception.

## Usage
To be written soon...

## Install
This is a single file, header-only library. Just place the ```ndarray.hpp```
file in your projects include directory, inorder to include it for use.
