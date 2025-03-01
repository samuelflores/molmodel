#ifndef SimTK_MOLMODEL_ATOM_H_
#define SimTK_MOLMODEL_ATOM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Christopher Bruns                                                 *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "molmodel/internal/Compound.h"
#include "molmodel/internal/AtomSubsystem.h"
#include <map>
#include <ostream>
namespace SimTK {


// Class to help manage indexing atoms, subcompounds, and bond centers
class CompoundPathName {
public:

    static bool isValidTokenName(const String& path) {
        // no empty strings allowed
        if (path.size() == 0) return false;

        for (int c = 0; c < (int)path.size(); ++c) {
            const char& ch = (const char&) path[c];

            // No "/" separator characters allowed in token name
            if ('/' == ch) return false;

            // Space characters are OK...
            // if ( ! isalnum(path[c]) ) return false;
        }

        return true;
    }

    static bool isValidAtomName(const String& path) {
        return isValidTokenName(path);        
    }

    static bool isValidSubcompoundName(const String& path) {
        return isValidTokenName(path);        
    }

    static bool isValidBondCenterName(const String& path) {
        return isValidTokenName(path);        
    }

    // Return the path with the left component removed
    static String shiftLeftPathName(const String& path) 
    {
        int start = (int)path.find_first_of("/");

        ++start; // step past "/" character
//std::cout << "CompundAtom path=" << path << "  start=" << start << std::endl;

        assert( start < (int)path.size() );

        return path.substr(start);
    }

    // Take pathname of the form xxx/yyy/zzz, check its validity and optionally
    // return as a list of separate subsystem names. We return true if we're successful,
    // false if the pathname is malformed in some way. In that case the last segment
    // returned will be the one that caused trouble.
    static bool isValidSubcompoundPathName(const String& pathname, 
                                           std::vector<String>* segments)
    {
        String t;
        const int end = (int)pathname.size();
        int nxt = 0;
        if (segments) segments->clear();
        bool foundAtLeastOne = false;
        // for each segment
        while (nxt < end) {
            // for each character of a segment
            while (nxt < end && pathname[nxt] != '/')
                t += pathname[nxt++];
            foundAtLeastOne = true;
            if (segments) segments->push_back(t);
            if (!isValidTokenName(t))
                return false;
            t.clear();
            ++nxt; // skip '/' (or harmless extra increment at end)
        }
        return foundAtLeastOne;
    }


};


class Bond {
public:
    Bond() :
            // amRotatable(true),
            mobility(BondMobility::Default),
            amRingClosingBond(false)
    {}

    Bond(mdunits::Length l, Angle t, bool isRingClosing) :
            defaultLength(l), defaultDihedral(t), 
            // amRotatable(true), 
            mobility(BondMobility::Default),
            amRingClosingBond(isRingClosing)
    {}

    bool isRingClosingBond() const {
        return amRingClosingBond;
    }

    Bond& setRingClosingBond(bool b) {
        amRingClosingBond = b;
        return *this;
    }


    BondMobility::Mobility getMobility() const {return mobility;}
    Bond& setMobility(BondMobility::Mobility m) 
    {
        mobility = m;
        return *this;
    }
    //bool isRotatable() const {
    //    return amRotatable;
    //}

    // Bond& setRotatable(bool b) {
    //     amRotatable = b;
    //     return *this;
    // }


    Bond& setDefaultDihedralAngle(Angle angle) {
        defaultDihedral = angle;
        assert( (defaultDihedral >= 0) || (defaultDihedral <= 0) );
        return *this;
    }

    Bond& setDihedralAngle(State& state, SimbodyMatterSubsystem& matter, Angle angleInRadians) {
        assert(pinJointId.isValid());

        MobilizedBody::Pin& pin = (MobilizedBody::Pin&) matter.updMobilizedBody(pinJointId);
        pin.setAngle(state, angleInRadians);

        return *this;
    }

    Angle getDefaultDihedralAngle() const {
        assert( (defaultDihedral >= 0) || (defaultDihedral <= 0) );
        return defaultDihedral;
    }

    Bond& setPinBody(MobilizedBody::Pin& pin) 
    {
        pinJointId = pin.getMobilizedBodyIndex();

        pin.setDefaultAngle(defaultDihedral);

        return *this;
    }

    Bond& setRiboseBody(MobilizedBody::FunctionBased& pin) 
    {
        pinJointId = pin.getMobilizedBodyIndex();

        // pin.setDefaultAngle(defaultDihedral);// TODO

        return *this;
    }

    MobilizedBodyIndex getPinJointId() const {return pinJointId;}

    Angle getDihedralAngle(const State& state, const SimbodyMatterSubsystem& matter) const {
        assert(pinJointId.isValid());

        const MobilizedBody::Pin& pin = (const MobilizedBody::Pin&) matter.getMobilizedBody(pinJointId);

        return pin.getAngle(state);
    }

    Bond& setDefaultBondLength(mdunits::Length d) {
        defaultLength = d;
        return *this;
    }

    mdunits::Length getDefaultBondLength() const {
        return defaultLength;
    }

    // This transform is symmetric, so it does not matter which bondCenterIsRelative to which
    Transform getDefaultBondCenterFrameInOtherBondCenterFrame() const {
        //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
        // 1) rotate about x-axis by dihedral angle
        Transform dihedral(Rotation(getDefaultDihedralAngle(), XAxis));

        // 2) translate along x-axis by bond length
        Transform bondLength(Vec3(defaultLength, 0, 0));

        // 3) rotate 180 degrees about y-axis to face the parent bond center
        Transform aboutFace( Rotation(180*Deg2Rad, YAxis) );

        Transform BC1_X_BC2 = dihedral * bondLength * aboutFace;
        //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;

        return BC1_X_BC2;
    }

	Angle getDefaultDihedral() const {return defaultDihedral;}

private:
    Angle defaultDihedral;
    mdunits::Length defaultLength;

    // bool amRotatable;
    BondMobility::Mobility mobility;

    bool amRingClosingBond;

    MobilizedBodyIndex pinJointId;
};


// BondInfo is Bond as viewed from a particular compound.
// These are indexed by Compound::BondIndex.
class BondInfo {
public:
    BondInfo() {}

    // 1 Ring closing bond - bonds that exist in addition to the tree structure of bonds
    BondInfo(
        Compound::BondIndex,
        Compound::BondCenterIndex bondCenter1,
        Compound::BondCenterIndex bondCenter2,
        const Bond& bond
        );

    bool isRingClosingBond() const {
        return bond.isRingClosingBond();
    }

    //bool isLocalBond() const {
    //    return amLocalBond;
    //}

    const Bond& getBond() const {
        // assert (isLocalBond() || isRingClosingBond());
        return bond;
    }

    Bond& updBond() {
        // assert (isLocalBond() || isRingClosingBond());
        return bond;
    }

    Compound::BondIndex getIndex() const {return id;}

    Compound::BondCenterIndex getParentBondCenterIndex() const {
        return parentBondCenterIndex;
    }
    Compound::BondCenterIndex getChildBondCenterIndex() const {
        return childBondCenterIndex;
    }

    BondInfo& setParentBondCenterIndex(Compound::BondCenterIndex i) {
        parentBondCenterIndex = i;
        return *this;
    }
    BondInfo& setChildBondCenterIndex(Compound::BondCenterIndex i) {
        childBondCenterIndex = i;
        return *this;
    }

private:
    Compound::BondIndex id;
    Compound::BondCenterIndex childBondCenterIndex;
    Compound::BondCenterIndex parentBondCenterIndex;

    // 3 possibilities of bond type
    // 1) Ring closing bond, originally defined in this compound (as opposed to within a subcompound)
    // bool amRingClosingBond; // the bond should know

    // 2) Subcompound bond - bond created within a subcompound - so delegate to subcompound
    //bool             amSubcompoundBond;
    //Compound::Index     subcompoundId;
    //Compound::BondIndex subcompoundBondIndex;
    
    // 3) Local bond - primary bond connecting this compound to an immediate subcompound
    // bool             amLocalBond;
    Bond             bond;
    //Compound         localBondSubcompound;

    // Compound childCompound;
    friend std::ostream& operator<<(std::ostream&, const BondInfo&);
};

// BondCenter is a class representing one half of a covalent
// bond between two atoms.  BondCenter belongs to one atom and
// represents one possible direction for a covalent bond.
class BondCenter {
public:
    enum Chirality {RightHanded, LeftHanded, Planar};

    // Use Paul's method
    // Given two other bond directions, and two bond angles,
    // determine direction of third bond
    static UnitVec3 getBondDirection(
        const UnitVec3& a1,
        Angle theta1,
        const UnitVec3& a2,
        Angle theta2,
        Chirality chirality
        ) 
    {
        //std::cout//<<__FILE__<<":"<<__LINE__<<" a1 theta1 a2 theta2 chirality "<< a1<<", "<< theta1<<", "<<a2<<", "<<theta2<<", "<<chirality<<std::endl;
        // Complete basis with third vector
        const UnitVec3 a3(a1 % a2);
        //std::cout//<<__FILE__<<":"<<__LINE__<<std::endl;

        if (chirality == Planar) {
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            assert( theta1 >= 0 );
            // assert( theta1 >= -180*Deg2Rad );
            assert( theta1 <= 180*Deg2Rad );

            assert( theta2 >= 0 );
            assert( theta2 <= 180*Deg2Rad );

            // TODO - decide direction of theta1 based on estimate of theta2

            // estimate angle between first two bond centers
            Angle theta12 = dot(a1, a2);
            if (theta12 > 1.0) theta12 = 1.0;
            if (theta12 < -1.0) theta12 = -1.0;
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            theta12 = std::acos(theta12);

            // is this bond on opposite side from bond2 with respect to bond1?
            // two cases, 1) 1-3 + 1-2 == 2-3
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            Angle oppositeError1 = std::abs(theta2 - theta1 - theta12);
            // 2) 1-3 + 1-2 == 360 - 2-3
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            Angle oppositeError2 = std::abs(360*Deg2Rad - theta12 - theta1 - theta2);

            // or is this bond on the same side as bond2? => 2-3 angle == difference between 1-3 and 1-2
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            Angle adjacentError = std::abs(std::abs(theta12 - theta1) - theta2);

            // if this is on the opposite side from bond2, flip theta1
            if ( (adjacentError > oppositeError1) || (adjacentError > oppositeError2) ) 
                theta1 = -theta1; // move to opposite hemisphere from a2 w.r.t. a1

            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            UnitVec3 direction( Rotation( theta1, a3 ) * a1 );
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            UnitVec3 sanityCheck( Rotation( theta12, a3) * a1 );

            return direction;
        }

        else { // non-planar chirality
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;

            // Bond angles are strictly positive
            assert(theta1 >= 0);
            assert(theta1 <= 180*Deg2Rad);
            assert(theta2 >= 0);
            assert(theta2 <= 180*Deg2Rad);


            // Compute coefficients of new bond direction in basis a1,a2,a3
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            Real cosTheta = dot(a1, a2);
            //std::cout<<__FILE__<<":"<<__LINE__<<" cosTheta "<<cosTheta<<std::endl;
            assert(cosTheta != 0); // non colinear
            assert(cosTheta >= -1);
            assert(cosTheta <= 1);

            Real cosTheta1 = std::cos(theta1);
            //std::cout<<__FILE__<<":"<<__LINE__<<" cosTheta1 "<<cosTheta1<<std::endl;
            Real cosTheta2 = std::cos(theta2);
            //std::cout<<__FILE__<<":"<<__LINE__<<" cosTheta2 "<<cosTheta2<<std::endl;
            Real sinSquaredTheta = 1.0 - (cosTheta*cosTheta);
            //std::cout<<__FILE__<<":"<<__LINE__<<" sinSquaredTheta "<<sinSquaredTheta<<std::endl;

            // a1 and a2 must not be parallel
            assert(sinSquaredTheta > 0);

            Real v1 = (cosTheta1 - cosTheta*cosTheta2) / sinSquaredTheta;
            Real v2 = (cosTheta2 - cosTheta*cosTheta1) / sinSquaredTheta;
            //std::cout<<__FILE__<<":"<<__LINE__<<" v1 v2 = "<<v1<<", "<<v2<<std::endl;
            Real v3Squared = 1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta);

	    //std::cout<<__FILE__<<":"<<__LINE__<<" v1 v2 cosTheta = "<<v1<<", "<<v2<<", "<<cosTheta<<std::endl; 
	    //std::cout<<__FILE__<<":"<<__LINE__<<" 1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta) = "<<1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta)<<std::endl;
            if ( chirality == LeftHanded )     {//std::cout<<__FILE__<<":"<<__LINE__<<" chiralit =  LeftHanded"<<std::endl; 
	    }
	    else if (chirality == RightHanded) {//std::cout<<__FILE__<<":"<<__LINE__<<" chiralit = RightHanded"<<std::endl; 
	    }
	    else                               {std::cout<<__FILE__<<":"<<__LINE__<<" unknown chirality     "<<std::endl; }
            // no solutions for certain sets of angles
            //assert(v3Squared >= 0);
            if (!(v3Squared >= 0)){
                std::cout<<__FILE__<<":"<<__LINE__<<" No solution found .. v1 v2 cosTheta = "<<v1<<", "<<v2<<", "<<cosTheta<<". Consider a larger planarity threshold."<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" 1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta) = "<<1.0 - (v1*v1 + v2*v2 + 2.0*v1*v2*cosTheta)<<std::endl;
		// If you want to read in spiral DNA, comment out this:
                exit(1);
		// And uncomment this:                                      
		//v3Squared = -v3Squared;
            }

            Real v3 = std::sqrt(v3Squared);

            if ( chirality == LeftHanded ) v3 = -v3;

            return UnitVec3(v1*a1 + v2*a2 + v3*a3);
        }
    }

    // Default constructor to be used only for storage in std types
    BondCenter();

    BondCenter(Angle angle1, Angle angle2, int yCenter, Chirality c);

    bool isBonded() const;

    BondCenter& setDefaultBondLength(mdunits::Length l) {
        assert (! isBonded() );

        defaultBondLength = l;

        return *this;
    }

    BondCenter& setDefaultDihedralAngle(Angle a) {
        assert (! isBonded() );

        defaultDihedralAngle = a;

        return *this;
    }

    mdunits::Length getDefaultBondLength() const {
        return defaultBondLength;
    }

    Angle getDefaultDihedralAngle() const {
        return defaultDihedralAngle;
    }

    int getDefaultDihedralReferenceCenter() const {
        return defaultDihedralReferenceCenter;
    }

    Angle getDefaultBond1Angle() const {return defaultBond1Angle;}
    BondCenter& setDefaultBond1Angle(Angle angle) 
    {
        ////std::cout<<__FILE__<<":"<<__LINE__  << " Problem is here! At the matchDefaultBondAngles stage"<<std::endl;
        assert(angle > 0);
        assert(angle <= SimTK::Pi);
        defaultBond1Angle = angle;
        //std::cout<<__FILE__<<":"<<__LINE__  << " angle "<< angle <<std::endl;
        return *this;
    }

    Angle getDefaultBond2Angle() const {return defaultBond2Angle;}
    BondCenter& setDefaultBond2Angle(Angle angle) 
    {
        assert(angle > 0);
        assert(angle <= SimTK::Pi);
        defaultBond2Angle = angle;
        return *this;
    }

    Chirality getChirality() const {return chirality;}

    BondCenter& setChirality(Chirality c) {
        chirality = c;
        return *this;
    }

    BondCenter& setInboard(bool b) {
        assert(!bonded);
        inboard = b;
        return *this;
    }

    bool isInboard() const {
        return inboard;
    }

    BondCenter& setBonded(bool b);

protected:

private:

    bool inboard;
    bool bonded;

    Angle defaultBond1Angle; // bond angle with first bond center (on x-axis)
    Angle defaultBond2Angle; // bond angle with second bond center (in +y half of XY plane)
    int defaultDihedralReferenceCenter; // index of other bond center to use as y-axis
    Chirality chirality;

    // Compound::BondIndex bondIndex; // if bonded, the bond object containing a subcompound

    mdunits::Length defaultBondLength;
    Angle defaultDihedralAngle;
};


class CompoundAtom {
public:
    // Atom::BondCenterIndex is a unique int index type. Don't confuse
    // it with Compound::BondCenterIndex -- this is a different index restricted to
    // a much smaller number.

    SimTK_DEFINE_UNIQUE_LOCAL_INDEX_TYPE(CompoundAtom,BondCenterIndex)
    
    CompoundAtom(const Element * element, Transform location = Transform());

    const Element * getElement() const;

    bool hasBondCenter(BondCenterIndex index) const {
        return 0 <= index && index < (int)bondCenters.size();
    }

    BondCenterIndex addFirstBondCenter() {
        assert( ! hasBondCenter(BondCenterIndex(0)) );

        // first bond center uses second bond center as dihedral reference
        int refY = 1;
        const BondCenterIndex bondCenterIndex = BondCenterIndex(bondCenters.size());

        bondCenters.push_back(BondCenter(0, 0, refY, BondCenter::Planar));

        assert( bondCenterIndex == 0 );
        assert( hasBondCenter(BondCenterIndex(0)) );

        return bondCenterIndex;
    }

    BondCenterIndex addSecondBondCenter(
        Angle bond1Angle
        ) 
    {
        assert( hasBondCenter(BondCenterIndex(0)) );
        assert( ! hasBondCenter(BondCenterIndex(1)) );

        assert( bond1Angle > 0 );
        assert( bond1Angle <= 180*Deg2Rad );

        // all bond centers except first use first bond center as dihedral reference
        int refY = 0;
        const BondCenterIndex bondCenterIndex = BondCenterIndex(bondCenters.size());

        bondCenters.push_back(BondCenter(bond1Angle, 0, refY, BondCenter::Planar));

        // TODO - ensure that theta1 remains in sync with second bond center angle
        // theta1 = bond1Angle;

        assert( bondCenterIndex == 1 );
        assert( hasBondCenter(BondCenterIndex(1)) );

        return bondCenterIndex;
    }

    BondCenterIndex addPlanarBondCenter(
        Angle bond1Angle,
        Angle bond2Angle
        ) 
    {
        assert( hasBondCenter(BondCenterIndex(0)) );
        assert( hasBondCenter(BondCenterIndex(1)) );

        assert( bond1Angle > 0 );
        assert( bond1Angle <= 180*Deg2Rad );

        assert( bond2Angle > 0 );
        assert( bond2Angle <= 180*Deg2Rad );

        // all bond centers except first use first bond center as dihedral reference
        int refY = 0;

        const BondCenterIndex bondCenterIndex = BondCenterIndex(bondCenters.size());
        assert( ! hasBondCenter(bondCenterIndex) );

        bondCenters.push_back(BondCenter(bond1Angle, bond2Angle, refY, BondCenter::Planar));

        assert( hasBondCenter(bondCenterIndex) );

        return bondCenterIndex;
    }

    BondCenterIndex addRightHandedBondCenter(
        Angle bond1Angle,
        Angle bond2Angle
        ) 
    {
        assert( hasBondCenter(BondCenterIndex(0)) );
        assert( hasBondCenter(BondCenterIndex(1)) );

        assert( bond1Angle > 0 );
        assert( bond1Angle <= 180*Deg2Rad );

        assert( bond2Angle > 0 );
        assert( bond2Angle <= 180*Deg2Rad );

        // all bond centers except first use first bond center as dihedral reference
        int refY = 0;

        const BondCenterIndex bondCenterIndex = BondCenterIndex(bondCenters.size());
        assert( ! hasBondCenter(bondCenterIndex) );

        bondCenters.push_back(BondCenter(bond1Angle, bond2Angle, refY, BondCenter::RightHanded));

        assert( hasBondCenter(bondCenterIndex) );

        return bondCenterIndex;
    }

    BondCenterIndex addLeftHandedBondCenter(
        Angle bond1Angle,
        Angle bond2Angle
        ) 
    {
        assert( hasBondCenter(BondCenterIndex(0)) );
        assert( hasBondCenter(BondCenterIndex(1)) );

        assert( bond2Angle > 0 );
        assert( bond2Angle <= 180*Deg2Rad );

        // all bond centers except first use first bond center as dihedral reference
        int refY = 0;

        const BondCenterIndex bondCenterIndex = BondCenterIndex(bondCenters.size());
        assert( ! hasBondCenter(bondCenterIndex) );

        bondCenters.push_back(BondCenter(bond1Angle, bond2Angle, refY, BondCenter::LeftHanded));

        assert( hasBondCenter(bondCenterIndex) );

        return bondCenterIndex;
    }

    BondCenter& updBondCenter(BondCenterIndex index) {
        return bondCenters[index];
    }
    const BondCenter& getBondCenter(BondCenterIndex index) const {
        return bondCenters[index];
    }

    CompoundAtom& setBiotypeIndex(BiotypeIndex bt) {
        biotypeIx = bt;
        return *this;
    }

    BiotypeIndex getBiotypeIndex() const {return biotypeIx;}

    const Transform& getDefaultFrameInCompoundFrame() const {
        return localTransform;
    }

    CompoundAtom& setDefaultFrameInCompoundFrame(const Transform& transform) {
        localTransform = transform;
        return *this;
    }

    int getNumBonds() const {
        return (int)bondCenters.size();
    }

    /** 
      @return the bond angle, in radians, between bond 1 and bond 2.
      returns NaN if there are less than two bonds, or if the angle is undefined
    */
    Angle getDefaultBond12Angle() const {
        if (getNumBonds() < 2) return NaN;
        else return getBondCenter(BondCenterIndex(1)).getDefaultBond1Angle();
    }

    UnitVec3 getBondCenterDirectionInAtomFrame(BondCenterIndex index) const 
    {
        const BondCenter& bondCenter = getBondCenter(index);

        if (index == 0) return UnitVec3(1, 0, 0);
        else if (index == 1) {
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            Angle theta = getBondCenter(index).getDefaultBond1Angle();
            //std::cout<<__FILE__<<":"<<__LINE__<<" theta "<<theta<< std::endl;

            UnitVec3 xAxis(1,0,0);
            Rotation rotMat(theta, ZAxis);
            UnitVec3 direction(rotMat * xAxis);

            return direction;
        }
        else {
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            UnitVec3 a1 = getBondCenterDirectionInAtomFrame(BondCenterIndex(0));
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            UnitVec3 a2 = getBondCenterDirectionInAtomFrame(BondCenterIndex(1));
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            Angle theta1 = bondCenter.getDefaultBond1Angle();
            Angle theta2 = bondCenter.getDefaultBond2Angle();
            BondCenter::Chirality chirality = bondCenter.getChirality();
            //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
            //std::cout<<__FILE__<<":"<<__LINE__<<" a1 theta1 a2 theta2 chirality "<< a1<<", "<< theta1<<", "<<a2<<", "<<theta2<<", "<<chirality<<std::endl;
            return BondCenter::getBondDirection(a1, theta1, a2, theta2, chirality);
        }
    }

    Transform calcDefaultBondCenterFrameInAtomFrame(BondCenterIndex bondCenterIndex) const 
    {
        //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
        // use getBondCenter() method to manage sanity checks
        const BondCenter& bondCenter = getBondCenter(bondCenterIndex);

        UnitVec3 direction = getBondCenterDirectionInAtomFrame(bondCenterIndex);

        // Another bondCenter is used to define the bond center y-axis for dihedral angles
        // Reference direction for dihedral computation
        const BondCenterIndex yAxisIndex = BondCenterIndex(bondCenter.getDefaultDihedralReferenceCenter());
        UnitVec3 ydir(0, 1, 0); // default to actual y-axis for one-bond-center-atom case
        // Some atoms have only one bond center, for others...
        if ((yAxisIndex != bondCenterIndex) && (yAxisIndex < (int)bondCenters.size()))
            ydir = getBondCenterDirectionInAtomFrame(yAxisIndex);

        // This creates a Rotation whose X axis is in "direction", and whose Y axis
        // is (at least roughly) in direction ydir.
        return Transform(Rotation(direction, XAxis, ydir, YAxis));
    }


    BondCenterIndex getInboardBondCenterIndex() const
    {
        for (BondCenterIndex b(0); b < getNumBonds(); ++b)
        {
            const BondCenter& center = getBondCenter(b);
            if (! center.isInboard()) continue;
            return b;
        }
        return BondCenterIndex(); // invalid index
    }

    Transform calcDefaultFrameInInboardCenterFrame() const
    {
        return ~calcDefaultBondCenterFrameInAtomFrame(getInboardBondCenterIndex());
    }
    
    // RUNTIME
    void setAtomSubsystemAtomIndex(AtomSubsystem::AtomIndex ssIx) {
        atomSubsystemAtomIndex = ssIx;
    }
    
    AtomSubsystem::AtomIndex getAtomSubsystemAtomIndex() const {
        return atomSubsystemAtomIndex;
    }
    
    void setDuMMAtomIndex(DuMM::AtomIndex dummId) {
        dummAtomIndex = dummId;
    }
    DuMM::AtomIndex getDuMMAtomIndex() const {
        return dummAtomIndex;
    }
    void setDuMMPrimaryClusterIndex(DuMM::ClusterIndex dummId) {
        dummPrimaryClusterIndex = dummId;
    }
    DuMM::ClusterIndex getDuMMPrimaryClusterIndex() const {
        return dummPrimaryClusterIndex;
    }

    Vec3 getLocationInMobilizedBodyFrame() const 
    {
        return frameInMobilizedBodyFrame.p();
    }
    const Transform& getFrameInMobilizedBodyFrame() const 
    {
        return frameInMobilizedBodyFrame;
    }
    CompoundAtom& setFrameInMobilizedBodyFrame(Transform loc) {
        frameInMobilizedBodyFrame = loc;
        return *this;
    }

    CompoundAtom& setMobilizedBodyIndex(MobilizedBodyIndex bodyIx) {
        mobilizedBodyIndex = bodyIx;
        return *this;
    }

    MobilizedBodyIndex getMobilizedBodyIndex() const {
        return mobilizedBodyIndex;
    }

	size_t getNumBondCenters() const {
		return bondCenters.size();
	}

private:
    const Element *             element;
//    int                         formalCharge;
    BiotypeIndex                   biotypeIx;
    Transform                   localTransform; // relative to the parent Compound's frame
    std::vector<BondCenter>     bondCenters;

    // RUNTIME MEMBERS
    AtomSubsystem::AtomIndex atomSubsystemAtomIndex;

    DuMM::AtomIndex dummAtomIndex;
    DuMM::ClusterIndex dummPrimaryClusterIndex; // proxy for rigid body id before bodies are actually specified

    MobilizedBodyIndex mobilizedBodyIndex; // populate after bodies are specified
    Transform frameInMobilizedBodyFrame;

    // Angle theta1;  // angle between bond centers one and two
    friend std::ostream& operator<<(std::ostream&, const CompoundAtom&);
};


// BondCenterInfo represents the view of an existing BondCenter,
// from the viewpoint of a particular parent compound
class BondCenterInfo {
public:
    typedef std::pair<Compound::AtomIndex, CompoundAtom::BondCenterIndex> AtomKey;
    typedef std::map<AtomKey, Compound::BondCenterIndex>          AtomKeyMap;

    BondCenterInfo() {
    }

    BondCenterInfo(
        Compound::BondCenterIndex  i,
        Compound::AtomIndex        a,
        CompoundAtom::BondCenterIndex   c) 
        : 
        id(i), 
        atomId(a), 
        atomBondCenterIndex(c)
    {}

    // const Compound::BondCenterName& getName() const;

    // BondCenterInfo& setName(const Compound::BondCenterName& n) {
    //     name = n;
    //     return *this;
    // }

    Compound::BondCenterIndex getIndex() const {return id;}

    Compound::AtomIndex getAtomIndex() const {return atomId;}

    CompoundAtom::BondCenterIndex getAtomBondCenterIndex() const {return atomBondCenterIndex;}

    const Compound::BondIndex getBondIndex() const {
        assert( isBonded() );
        return bondIndex;
    }

    BondCenterInfo& setBondIndex(Compound::BondIndex b) {
        assert( ! isBonded() );
        bondIndex = b;
        assert(isBonded());
        return *this;
    }

    Compound::BondCenterIndex getParentCompoundBondCenterIndex() const {
        return parentCompoundBondCenterIndex;
    }

    BondCenterInfo& setParentCompoundBondCenterIndex(Compound::BondCenterIndex i) {
        parentCompoundBondCenterIndex = i;
        return *this;
    }

    Compound::BondCenterIndex getBondPartnerBondCenterIndex() const {
        return bondPartnerBondCenterIndex;
    }

    BondCenterInfo& setBondPartnerBondCenterIndex(Compound::BondCenterIndex i) {
        bondPartnerBondCenterIndex = i;
        return *this;
    }

    bool isBonded() const { 
        return bondIndex.isValid();
    }


private:
    Compound::BondCenterIndex id;
    Compound::BondCenterIndex bondPartnerBondCenterIndex;
    Compound::BondCenterIndex parentCompoundBondCenterIndex;

    Compound::BondIndex       bondIndex; // if bonded

    // These map this BondCenterInfo to the particular (AtomInfo,atom bond center#) pair to
    // which it refers.
    Compound::AtomIndex       atomId;
    CompoundAtom::BondCenterIndex  atomBondCenterIndex;
};

// These are indexed by Compound::AtomIndex.
class AtomInfo {
public:
    // One constructor for local atoms
    AtomInfo(Compound::AtomIndex id, const CompoundAtom& atom, bool isBase /*, const Transform& frameInCompound = Transform()*/ );

    AtomInfo(const AtomInfo &other)
        : id(other.id),
          name(other.name),
          synonyms(other.synonyms),
          bIsBaseAtom(other.bIsBaseAtom),
          atom(other.atom)
    {}

    AtomInfo(AtomInfo &&other) noexcept
        : id(std::move(other.id)),
          name(std::move(other.name)),
          synonyms(std::move(other.synonyms)),
          bIsBaseAtom(other.bIsBaseAtom),
          atom(std::move(other.atom))
    {}

    bool hasValidAtomName() const {
        return CompoundPathName::isValidAtomName(name);
    }

    const Compound::AtomName& getName() const {
        // assert(CompoundPathName::isValidAtomName(name));
        return name;
    }

    const std::set<Compound::AtomName>& getNames() const {
        return synonyms;
    }

    const AtomInfo& addName(Compound::AtomName n) {
        // assert(CompoundPathName::isValidAtomName(n));
        name = std::move(n);
        synonyms.insert(name);
        return *this;
    }

    bool isBaseAtom() const {return bIsBaseAtom;}

    Compound::AtomIndex getIndex() const {return id;}

    const CompoundAtom& getAtom() const {return atom;}
    CompoundAtom& updAtom() {return atom;}

private:
    Compound::AtomIndex id;

    // Compound::AtomIndex parentCompoundAtomIndex;

    Compound::AtomName name;
    std::set<Compound::AtomName> synonyms;

    bool bIsBaseAtom;
    // Transform frameInCompound;

    CompoundAtom atom;

    friend std::ostream& operator<<(std::ostream&, const AtomInfo&);
};


} // namespace SimTK

#endif // SimTK_MOLMODEL_ATOM_H_
