# Rigid body simulation and control journey


## Summary

- Writing a simulation engine from scratch Python
- Reusing the engine for solving optimal control tasks
- Enjoying the results in Blender



## Story of none other than my life

In the usual depth first search way:
- Need to learn about attitude control - I like rockets
- Hmm, multiple reaction wheels with two motors each, complicated
- Oops, I don't have a good understanding of a simple gyroscope (spinnig wheel -_-)
- Let's make a general purpose rigid body simulation engine that can handle many type of joints!

To make progress on these problems, few constraints:
- I don't want to write equations for each new problem I study
- matlab is out - never managed to reuse anything I've written in it (and simulink :gasp:)
- python is in:
    - sunk cost fallacy - big chunk of my personal library is written in it
    - the ecosystem is just too good
      - scientific libraries: numpy/scipy/pygmo
      - blender for visualization: with python scripting + jupyter integration it's a treat!
- I don't like doing the same thing quite, it'd be convenient to reuse the simulation engine for solving the optimal control tasks


There's a fairly strong incompatibility here: python - simulation engine - optimal control. Its name is speed.
After some failed attempt to use Sympy to solve this issue, I stumbled upon the awesome Jax project that solved this problem exactly how I wanted it.



## Rigid body engine


I'll be rather terse here, go read Featherstone's Rigid body simulation algorithms, first 6 chapters if you want a detailed explanation of how the code works.

Please keep in mind that I don't know much on the subject:
- I haven't read any other book
- went for the more familiar newtonian mechanics


Still, I found the notion of spatial vector presented in this book to be quite golden.

Writing newton laws for a 3D rigid body:
- F = dp/dt: F force, p is linear momentum
- T = dL/dt: T torque, L is angular momentum

The properties of a 3D rigid body can mostly be grouped into pairs, linear and angular:
- position: center of mass + orientation
- linear + angular speed/accelation
- linear momentum, angular momentum


Featherstone introduces objects with very nice properties that covers both the linear and angular components of the 3D rigid body:
- spatial vector for linear+angular speed
- spatial acceleration for linear+angular acceleration
- spatial force for force+torque
- spatial momentum for linear+angular momentum
- full inertia matrix for mass+inertia tensor


Time derivative, addition, change of basis are operations that needs to be defined for these objects, details are in the book.
Some equations we end up with:
- F = Ia + [v] Iv
    - F,a,v spatial force, acc, velocity. I spatial momentum
    - [v] is a differention operator because the spatial force Iv is rotating
- v_i = X_i_rw v_prev + qd: 
    - speed of rigid body i defined from the rigid body it is previously linked to (with X_i_rw the transform matrix from the previous body to the joint)
    - qd the speed of the joint (ie relative movement between the two linked bodies)
- a_i = X_i_rw a_prev + qdd + [vi] qd


It took quite a while to assimilate imperfectly some of this material (writing about it, I find some of my understanding have faded).


With these equations, we obtain an inverse dynamic primitive (force = f(acceleration)) which we invert to get a forward dynamics (acc=f(force)). Integrating, we have our simulation algorithm.


### Boring stuff

- Momentum conservation: without external forces, the integration scheme will cause the momentum of the system to slowly drift. Some additionnal mechanism is required to guarantee momentum conservation.
- jax
    - jax arrays can almost be used as a drop-in replacement of numpy arrays. The major source of adjustements comes from the immutability of the jax arrays. The jax docs goes into detail about it.
    - I want to use jax arrays only when necessary - no heavy-handed sed -> need a dispatcher for some operations


## Control theory


Once the engine is working fast enough, it's possible without many modifications to experiment with the controls of a model to have it reach a given state.


Using a general purpose optimization solver as a black box (not trying to model the problems as convex or whatever), we can inefficiently solve control problems as such:
- min_u f(u,start_state, end_state), u is the input forces we can provide on the system over time
- u = (u_1, ..., u_n) the input forces at time (t_i)
- s_(i+1) = engine(s_i, u_i), with s_1 = start_state

For instance, we could have f = d(end_state, s_(n+1)) for some distance function.


And basically that's it. All that's left is to choose which solver to use.
Considering python's speed it's not suprising that it is not the most popular language for trajectory optimization libraries.

Still I ended up finding ESA's pygmo library that wraps many existing solvers, don't know if it's a bad choice or no but sometimes it gives nice results (considering the no-effort modelling of the optimization task).


## Results

Here will be presented a few simulation/control problems to put all the tools developed to work.
Rendering is done by blender, controlled through it's python API.

Obviously, there never was any meta-parameters tuning necessary to have half-decent results ;)


### Gyroscope effect

- spherical joint for string-axis
- pivot for axis-wheel
- gravity

### T handle in space: Dzhanibekov effect

- free joint
- no gravity


### Moving a box from A to B

- 1 degree of freedom, translation along given axis for box-world

### Inverse pendulum

- translation joint for cart-world
- pivot joint for pendulum-cart
- gravity
- control force only on translation joint, null on the pivot



<details>
<summary>Scene buildup code</summary>
```python
  sctx = SceneContext()
  tx = RBTree(sctx=sctx, split_rigid=0)
  root = tx.add(
      RBDescEntry(
          data=RBData(base_name='root'),
          spec=SolidSpec.Box(3, 1, 1, 3),
          link_data=LinkData(
              spec=LinkSpec(
                  type=RigidBodyLinkType.XLT_Z,
                  wr=Transform.From(rot=make_rot_ab(x=Vec3.Z(), y=Vec3.Y()))
              )
          ),
      )
  )

  axis_l = 10
  axis = tx.add(
      RBDescEntry(
          data=RBData(base_name='axis'),
          link_data=LinkData(
              spec=LinkSpec(
                  type=RigidBodyLinkType.PIVOT_Z,
                  wr=Transform.From(rot=make_rot_ab(y=Vec3.Z(), x=Vec3.Y()), pos=[0, 0, 0]),
                  rl=Transform.From(pos=[0, -axis_l / 2, 0])
              )
          ),
          spec=SolidSpec.Box(0.2, 0.5, axis_l, 0.2),
          parent=root,
      )
  )

  tx.add(
      RBDescEntry(
          data=RBData(),
          link_data=LinkData(
              spec=LinkSpec(
                  type=RigidBodyLinkType.RIGID, wr=Transform.From(pos=[0, -axis_l / 2, 0])
              )
          ),
          spec=SolidSpec.Sphere(2, 1),
          parent=axis,
      )
  )
  res = tx.create(root)

  def ctrl2model(sim, ctrl):
    res = np.zeros(sctx.sys_spec.ctrl_packer.pos)
    res = g_oph.set(res, ctrl[0:1])[0:1]  # z xlt
    return res

  fm = ForceModel(
      nctrl_f=lambda n_: 1,
      model2ctrl=lambda sim, model: model[0:1],
      ctrl2model=ctrl2model,
      bias_force_f=bias_force_g,
  )
  return SceneData(sctx=sctx, fm=fm)

```
</details>

### Attitude control with reaction wheels

- free joint for body-world
- 3 pivot joint for each wheel-body, along different axes
- control force only on pivot joints, null on the free joint




