Diffractometer
==============

This class provides the transformation from the setting angles of a six circle diffractometer into reciprocla space. There for we use mainly four frames in reciprocal Space: :math:`\theta` -frame, :math:`\phi` -frame, :math:`U` -frame and :math:`(H,K,L)` -frame.

Lab Frame
---------

We define the lab frame :math:`\Sigma` like in the paper about the four circle diffraction setup (Busing1967), i.e. the :math:`Z` -axis up, the :math:`Y` -axis along the X-ray beam and the :math:`Y` -axis forms a right handed orthonganl system with the others.

Rotations and Setting Angles
----------------------------

The rotation by :math:`\mu` along the :math:`+Z`-axis leads to :math:`\Sigma'` followed by the angles familiar from the four circle setup, i.e.
rotation by :math:`\theta` along the :math:`+X'` -axis leads to :math:`\Sigma''` 
(:math:`\theta` -frame),
rotation by :math:`\chi`   along the :math:`+Y''` -axis leads to :math:`\Sigma'''`, and
rotation by :math:`\phi`   along the :math:`+X'''` -axis leads to :math:`\Sigma''''`
(:math:`\phi` -frame)
The detector rotations start in :math:`\Sigma'` with 
rotation by :math:`\delta` along the :math:`+X'` -axis into :math:`\Sigma^*`, and
rotation by :math:`\gamma` along the :math:`+Z^*` -axis into :math:`\Sigma^{**}`.


Diffractometer Class
--------------------

.. automodule:: pyspec.diffractometer
   :members:
