#ifndef __OBJECT_H__
#define __OBJECT_H__

#include <vector>

// forward declarations
namespace Fluid
{
	class Object;
	class Sphere;
	class Box;
	class Group;
}

#include "config.h"
#include "vector3d.h"
#include "fluid.h"

namespace Fluid
{

/** Objects that interact with the fluid simulation.
 * All objects that interact with the fluid simulation are inherited from 
 * this base class. They all have in common that they can be sources,
 * forces and solid obstacles. Perhaps this class will provide grid
 * resolution independent simulation settings one day, so be careful
 * with your coordinate system: Do not depend on the fact that object
 * space equals simulation space at this moment.
 *
 * \ingroup core
 *
 * \version $Revision: 1.6 $
 *
 * \date $Date: 2003/02/07 10:54:56 $
 *
 * \author Malte
 *
 * \todo 
 *
 * \bug 
 *
 */
	class Object
	{
	public:
	/** Initializes to object to be neutral (no source, no force, not solid)
	 */
		Object();
		virtual ~ Object()
		{														/* do nothing */
		};

	/** Sets the force vector of the object.
	 * 
	 * \param force force vector (direction and strength)
	 */
		void SetForce(const Vector3D force);

	/** Sets the source rate of the object.
	 * 
	 * \param source source rate (new particles/time)
	 */
		void SetSource(const float source);

	/** Specifies wether object is solid or not.
	 * Solid objects block the flow: Particles cannot diffuse into solid 
	 * objects and cannot flow through them. This also applies to the velocity.
	 * 
	 * \param solid solid state
	 */
		void SetSolid(const bool solid);

	/** Initializes the object according to the solver settings.
	 * Initialize size, sources, forces etc. to match the grid coordinates of
	 * the fluid solver. This step is required everytime you change one of
	 * these parameters.
	 *
	 * \todo check wether the calculations for Source and Force are correct.
	 * They seem to be, but I'm still unsure.
	 *
	 * \param solver The solver in which the object will be used
	 */
		virtual void Init(const Solver & solver);

	/** Tests wether a vertex lies inside the object.
	 * 
	 * \param pos vertex position
	 * \return true if inside the object
	 */
		virtual bool Inside(const VecInt3D pos) const = 0;

	/** Tests wether a vertex lies inside a solid object.
	 * 
	 * \param pos vertex position
	 * \return true if inside a solid object
	 */
		virtual bool InsideSolid(const VecInt3D pos) const;

	/** Retrieve the force at a given position.
	 * 
	 * \param pos vertex position
	 * \return force vector
	 */
		virtual const Vector3D GetForce(const VecInt3D pos) const;

	/** Retrieve the source rate at a given position.
	 * 
	 * \param pos vertex position
	 * \return source rate
	 */
		virtual float GetSource(const VecInt3D pos) const;

	/** Set change state.
	 * Used to reset change state or to explicitly marking changes (should 
	 * never be required).
	 * 
	 * \param changed change state
	 */
		virtual void SetChanged(bool changed);

	/** Check wether object changed since the last call to SetChanged( false ).
	 * 
	 * \return true if changed
	 */
		virtual bool GetChanged() const;


		virtual void Translate(float x, float y, float z){} // changed by tobias (23.01.03)

		virtual const char * getName() { return "Object";}

	private:
		 Vector3D mForce;						///< force vector (direction and strength)
		float mSource;							///< source rate (new particels / time )

		Vector3D mForceGrid;
		float mSourceGrid;

		bool mSolid;								///< object is solid, see SetSolid()

		bool mChanged;							///< track changes
		
		const Solver* mpSolver;
		
		static int mMaxID;	// changed by tobias (23.01.03)
		int mID;						// changed by tobias (23.01.03)

	public:
		int GetID() {return mID;} // changed by tobias (23.01.03)
		bool GetSolid() const { return mSolid; }
	};

/** A sphere object.
 *
 * \ingroup core
 *
 * \version $Revision: 1.6 $
 *
 * \date $Date: 2003/02/07 10:54:56 $
 *
 * \author Malte
 *
 * \todo
 *
 * \bug
 *
 */
	class Sphere:public Object
	{
	public:
	/** Creates a sphere.
	 *
	 * \param center center point
	 * \param radius surface radius
	 */
		Sphere(const Vector3D center, const float radius);

		virtual bool Inside(const VecInt3D pos) const;

		virtual void Init(const Solver & solver);

		virtual void Translate(float x, float y, float z); // changed by tobias (23.01.03)
		
		virtual const char * getName() { return "Sphere";}
	private:
		 Vector3D mCenter;					///< center point
		float mRadius;							///< surface radius

		VecInt3D mCenterGrid;				///< center point in grid coordinates
		int mRadius2Grid;						///< surface radius in grid space, squared for faster access

	public:
		const Vector3D& GetCenter() const { return mCenter; }
		float GetRadius() const { return mRadius; }
	};

/** A box object.
 * Boxes are always aligned to the coordinate system axis.
 *
 * \ingroup core
 *
 * \version $Revision: 1.6 $
 *
 * \date $Date: 2003/02/07 10:54:56 $
 *
 * \author Malte
 *
 * \todo
 *
 * \bug
 *
 */
	class Box:public Object
	{
	public:
	/** Creates a box.
	 * All points with coordinates between lt and gt lie inside the box.
	 *
	 * \param lt lower limit
	 * \param gt upper limit
	 * \return
	 */
		Box(const Vector3D lt, const Vector3D gt);

		virtual bool Inside(const VecInt3D pos) const;

		virtual void Init(const Solver & solver);
		
		virtual void Translate(float x, float y, float z); // changed by tobias (23.01.03)
		
		virtual const char * getName() { return "Box";}

	private:
		 Vector3D mCornerLT;				///< lower limit
		Vector3D mCornerGT;					///< upper limit

		VecInt3D mCornerLTGrid;			///< lower limit in grid coordinates
		VecInt3D mCornerGTGrid;			///< upper limit in grid coordinates
	};

/** A object group.
 * Object groups are used to create object hierarchies. You can add several
 * child objects to a group. All object tests are performed on the children,
 * which can be groups themselves.
 *
 * \ingroup core
 *
 * \version $Revision: 1.6 $
 *
 * \date $Date: 2003/02/07 10:54:56 $
 *
 * \author Malte
 *
 * \todo 
 *
 * \bug 
 *
 */
	class Group:public Object
	{
	public:
		typedef std::vector < Object * >ObjList;

	/** Return iterator to the front of the list of contained elements.
	 */
		ObjList::iterator begin() { return mElements.begin(); }
	/** Return iterator to the end of the list of contained elements.
	 */
		ObjList::iterator end() { return mElements.end(); }

	/** Add another object to the group.
	 *
	 * \param object pointer to the object
	 */
		void Add(Object * object);
	
	/** Remove object from group.
	 * \param it iterator referencing to the object
	 */
		void Remove( ObjList::iterator it ) { mElements.erase( it ); }

		virtual bool Inside(const VecInt3D pos) const;

		virtual bool InsideSolid(const VecInt3D pos) const;

		virtual const Vector3D GetForce(const VecInt3D pos) const;

		virtual float GetSource(const VecInt3D pos) const;

		virtual bool GetChanged() const;

		virtual void Init(const Solver & solver);

	private:
		ObjList mElements;					///< list of child objects

	public:
		const ObjList* ExposeObjectList() const { return &mElements; } // dirty hack (sorry)
	};

};															// namespace Fluid


#endif													// __OBJECT_H__