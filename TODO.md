##TODO

- [ ] Feature match asessment/improvement.  Let C be the 3d point (a part of the model) with its photometric descriptors (SIFT, etc) and A, B be the features in the 1st and 2nd images respectively:
```
       C
       /\
      /  \
     /    \
    /      \
   A<------>B
```

Current match procedure works like this  
  * A vs. B are matched based on the L2 (or other) distance with epipolar constraint
  * A and B independently matched back to the model and the match is considered a match only if both matched the same 3d point C.

We discussed 2 ways to improve this matching  
  * If any one of the links is missing it is possible to complete it based on the available information.  I will try to do this and analyze if this gives us a better match.
  * Filter the matches based on the reprojection error (as a result of rigid motion that is estimated in the end of the algorithm).  Again, I will analyze the performance.


- [ ] Discovery of new features:
* For each stereo pair keep along the list of 3d points that did not match anything in the dictionary
* Given a new pair of images (and the 3d of its feature points) we may verify which 3d points comply with the calculated 3d motion.  If we saw some points comply in a number of frames, we may add it to the model.

- [ ] Normal calculation. Normals may be used to uderstand whether this current feature points is currently visible or it is not.  Approx the shape locally by a plane that passes through current point and some of its neighbors.

6/10/2014

We spoke mostly about error modeling in 3d estimation.  One paper that covers the things we spoke about is this: [Error Modeling in Stereo Navigation](http://repository.cmu.edu/cgi/viewcontent.cgi?article=2592&context=compsci). Just to recap basic points of the paper:
- Image measurements are noisy and thus depth estimation has errors as well
- Previous work shows that modeling depth error as a scalar leads to poor results.  This happens because the uncertainty introduced by triangulation is not a simple scalar function of distance to the point; it is also skewed and oriented.  Nearby points have a fairly compact uncertainty while distant points have a more elongated uncertainty that is roughly aligned with the line of sight to the point
- Error distribution is modeled by 3d Normal distribution. Its covariance matrix is given by: V=Jdiag(Vl,Vr)J', where Vl,Vr are covariances of the image coordinates, which are unknown, but are assumed to be uncorrelated with variances of 1 (or few) pixels.  So, the error is modeled as a co-variance matrix, as we spoke.
- The authors also propose a way to solve rigid motion with the above error modeling in mind.  They derive a closed form solution for transnational motion and an iterative solution for motion that also have rotations (basically they use the estimated co-variance matrix as a norm to weight errors). One thing is that their way does 3d-to-3d, while, I think, most works today tend to do direct error minimization, so it would be nice to find some way to incorporate these error measurements into re-projection minimization (BTW if you know of some work that does it I'd be happy to get a reference).
- The authors also propose a way to incorporate this error modeling in the path estimation (Kalman Filter based, by I did not get into the details here, since I'm not close to that yet)
