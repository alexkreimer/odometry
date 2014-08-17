pose3d
======

Estimate pose of articulated objects from images

General description
-------------------

We create **3d models augmented with photometric information** which may be used for object detection and/or pose estimation. (Note: there is only basic articulation support in the algorithm just yet. Our emphasis is on 3d models for rigids).

By **model** we mean a collection of 3d object points with corresponding appearance information (e.g., SIFT descriptors for these points)


Below is an outline of the algorithm with concise explanation of each step.
- *Model initialization* is partly manual: the user marks rigid parts in each init frame (init frames are predefined).  These marks are used to init the model. We calculate 3d points (by triangulating the matches) and store their appearance. At the moment appearance is a SIFT descriptor, but other options are possible (for piecewise planar objects covariant descriptors may be a good option, need a ref here).  Dictionary is a collection of 3d object points with a list of appearance descriptors attached to each point.
- *Step* of pose estimation/model refinement. Givean a new (calibrated) pair of images do:
  1. Calculate feature points (currently SIFT) in these input images (denote by A,B the set of features in the left and the right images respectively)
  2. Match features from A to features in B.  The match is performed by choosing a nearest neighbor based on L2 (or other) metric with an addition of epipolar geometry constraint.
  3. Match current frames' features back to the model (Lowe's 2nd best). Keep only features that close a circle (i.e. features from both of the current frames match the same feature in the dictionary).  This helps to remove wrong matches
  4. Reconstruct 3d for features that matched using the calibration of the rig
  5. Solve rigid motion between current vs. dictionary and thus obtain current pose​ of the object.  This uses RANSAC and returns R,t and an inlier set.
  6. *Model update:* add new SIFT descriptors to correspoding feature points in the dictionary
