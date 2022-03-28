# qpm
Quaternions for coordinate transformation and rotation in paleomagnetism

Python scripts are prepared to perform coordinate transformation and rotation based on quaternions. For transformation from a sample coordinate to the in situ geographic coordinate, python functions for block samples oriented by (strike, dip) ("qisb"), drill cores oriented by (azimuth, plunge) ("qisd"), and horizontally unoriented vertical piston cores ("qisp") are available. Functions for tilt correction by (strike, dip) of the formation ("qtilt") and for VGP calculation ("qvgp") are also included. After installing numpy-quaternion, these scripts can be executed in python by importing "qpm" and following the example codes given as comments.
