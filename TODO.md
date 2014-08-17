#TODO

##[ ] Feature match asessment/improvement.
  Let XX be the 3d point (a part of the model) with its photometric descriptors (SIFT, etc).  Let A, B be the features in the 1st and 2nd images respectively:
         XX
         /\
        /  \
       /    \
      /      \
     A<------>B
  Current match procedure works like this: A vs. B are matched based on the L2 (or other) distance with epipolar constraint.  After this both A and B independently matched back to the model and the match is considered a match only if both matched the same 3d point XX.
  We discussed 2 ways to improve this matching:
  1. If any one of the links is missing it is possible to complete it based on the available information.  I will try to do this and analyze if this gives us a better match.
  2. Filter the matches based on the reprojection error (as a result of rigid motion that is estimated in the end of the algorithm).  Again, I will analyze the performance.

##[ ] Discovery of new features:
  1. For each stereo pair keep along the list of 3d points that did not match anything in the dictionary
  2. Given a new pair of images (and the 3d of its feature points) we may verify which 3d points comply with the calculated 3d motion.  If we saw some points comply in a number of frames, we may add it to the model.
