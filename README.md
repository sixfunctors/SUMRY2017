# SUMRY2017
Code for the enumeration of Hassett Chambers and associated constants. 
Developed by Connor Halleck-Dube as part of SUMRY 2017, in collaboration with
Kenneth Ascher, Daniel Gershenson, and Elaine Hou. Adapted from code by Nicolle Gruzling.

This body of code can be used to enumerate Goldilocks linear threshold functions on n variables, 
or equivalently, chambers of the Hassett decomposition of the moduli space of weighted stable curves.
It was used to obtain numerical results on the number of such chambers for n <= 9, which is
presented in (https://arxiv.org/abs/1709.03663). 

This code is adapted from code for the enumeration of all linear threshold functions, 
written by Nicolle Gruzling. Specific information on code attribution can be found in the respective source files.

This code makes use of the C++ Big Integer library written by Matt McCutchen, which is in the public domain (https://mattmccutchen.net/bigint/). It also makes use 
of the "concurrentqueue" implementation of a locking multi-producer, multi-consumer, thread-safe queue. This can be found at (https://github.com/cameron314/concurrentqueue), and is published under Simplified BSD license. 
