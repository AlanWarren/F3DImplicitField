#include <iostream>
#include <math.h>
#include <string>

#include <fnmatch.h>

#include <boost/pointer_cast.hpp>
#include <boost/foreach.hpp>

#include <Field3D/DenseField.h>
#include <Field3D/SparseField.h>
#include <Field3D/InitIO.h>
#include <Field3D/Field3DFile.h>
#include <Field3D/FieldInterp.h>

#include "ImplicitField.h"
#include "rx.h"
#include "RixInterfaces.h"

#define DEBUG

using namespace std;
using namespace Field3D;

bool matchString(const string&, const vector<string>&);

template <typename T>
typename Field<T>::Ptr getField(Field3DInputFile&, vector<string>&, vector<string>&, vector<string>&, 
        vector<string>&, const string&, const string&, bool);

template <class T>
class VertexField : public ImplicitVertexValue
{
    public:
        VertexField(Field3DInputFile&, FieldMapping*, LinearFieldInterp<T>, 
                vector<string>&, M44d, const string&);
        ~VertexField() {}

        RtFloat EvalField(const RtPoint p);
        void GetVertexValue(RtFloat *result, const RtPoint p) 
        {
            result[0] = EvalField(p);
        }
    private:
        M44d t_oTw;
        const string t_field;
        vector<string> t_partitions;
        vector<string> t_fieldName;
        vector<string> t_attribName;
        vector<string> t_scalarLayers;
        FieldMapping *t_mapping;
        LinearFieldInterp<T> t_interpolator;
        typename Field<T>::Ptr t_buffer;

};

template <class T>
VertexField<T>::VertexField(Field3DInputFile &file, FieldMapping *mapping, LinearFieldInterp<T> interp,
        vector<string> &partitionsVec, M44d oTw, const string &fieldName) : 
    t_mapping(mapping), 
    t_interpolator(interp), 
    t_partitions(partitionsVec),
    t_oTw(oTw),
    t_field(fieldName)
{
    // XXX: Test for vector values
    typedef typename Field<T>::Ptr SFieldList;
    typedef typename Field<FIELD3D_VEC3_T<T> >::Vec VFieldList;

    SFieldList sFields = file.readScalarLayers<T>(t_field);
    VFieldList vFields = file.readVectorLayers<T>(t_field);
    bool isVec = false;

    if (vFields.size() > 0) {
        isVec = true;
        t_buffer = getField<T>(file, t_partitions, t_fieldName, t_attribName, 
                            t_scalarLayers, t_field, t_field, isVec);
    } else {
        t_buffer = getField<T>(file, t_partitions, t_fieldName, t_attribName, 
                            t_scalarLayers, t_field, t_field, isVec);
    }

}

template <class T>
RtFloat VertexField<T>::EvalField(const RtPoint p) 
{
    float ret;
    V3d wsP(p[0], p[1], p[2]);
    V3d vsP;

    t_mapping->worldToVoxel(wsP, vsP);
    // XXX: Add out-of-bounds checks for accurate sampling
    return t_interpolator.sample(*t_buffer, vsP);
}

class F3DImplicitField : public ImplicitField
{
    public:
    	F3DImplicitField(string path, float blur, float bbox_mod);
    	~F3DImplicitField(){}
    
    	virtual RtFloat Eval(const RtPoint p);
        virtual void EvalMultiple(int neval, float *result, int resultstride, const RtPoint *p);
        virtual void GradientEval(RtPoint result, const RtPoint p) {result[0] = 0; result[1] = 0; result[2] = 0;}
    	virtual void Range(RtInterval result, const RtPoint corners[8], const RtVolumeHandle h);
        virtual void Motion(RtPoint result, const RtPoint p);
        virtual void BoxMotion(RtBound result, const RtBound b);
        virtual void MotionMultiple(int neval, RtPoint *result, const RtPoint *p);
        ImplicitVertexValue *CreateVertexValue(const RtToken name, int nvalue);
        V3d xformPoint(V3d &Pt, FieldMapping *mapping, int frustrum, int fromobject);

        template <typename T>
        vector<typename Field<T>::Ptr> VelocitySamples(vector<string> &velFieldN, vector<string> &velAttrib, 
                                                       vector<string> &sLayers, bool isVec);
        template <typename T>
        void ReadLayersAndSetupFields(typename Field<T>::Vec sFields);

        template <typename T>
        VertexField<T>* DeployField(const RtToken &pname, int &nvalue);

        template <typename T>
        void collectVelocityBounds(typename Field<T>::Ptr xbuf,
                                   typename Field<T>::Ptr ybuf,
                                   typename Field<T>::Ptr zbuf);
        template <typename T>
        void collectVectorVelocityBounds(typename Field<T>::Ptr vec);

    private:
        vector<string> m_partitions;
        vector<string> m_fieldName;
        vector<string> m_scalarLayers;
        float m_xmax, m_ymax, m_zmax;
        float m_xmin, m_ymin, m_zmin;
        float m_dxmax, m_dymax, m_dzmax;
        Box3i m_extents;
        Box3i m_datawin;
        string m_dataType;
        float m_blur; 
        float m_bbox_mod;
        // space transformation matrix
        M44d Mp;
        // message display
	    RixMessages *msgs;
        // field3d
        Field3DInputFile m_in;
        FieldMapping *m_mapping;
        LinearFieldInterp<float> m_finterpolator;
        LinearFieldInterp<double> m_dinterpolator;
        LinearFieldInterp<half> m_hinterpolator;
        // density buffers
        Field<float>::Ptr m_fbuffer;
        Field<double>::Ptr m_dbuffer;
        Field<half>::Ptr m_hbuffer;
        // velocity buffers
        Field<float>::Ptr m_velx_fbuffer;
        Field<float>::Ptr m_vely_fbuffer;
        Field<float>::Ptr m_velz_fbuffer;
        Field<double>::Ptr m_velx_dbuffer;
        Field<double>::Ptr m_vely_dbuffer;
        Field<double>::Ptr m_velz_dbuffer;
        Field<half>::Ptr m_velx_hbuffer;
        Field<half>::Ptr m_vely_hbuffer;
        Field<half>::Ptr m_velz_hbuffer;
        // velocity mappings
        FieldMapping *m_velx_mapping;
        FieldMapping *m_vely_mapping;
        FieldMapping *m_velz_mapping;
        // bounds modifiers
        float m_vel_xmax, m_vel_ymax, m_vel_zmax;
        float m_vel_length;
}

F3DImplicitField::F3DImplicitField(string path, float blur, float bbox_mod) : 
    m_blur(blur), m_bbox_mod(bbox_mod)
{
    typedef Field3D::half half;

	RixContext *rixCtx = RxGetRixContext();
	msgs = (RixMessages*)rixCtx->GetRixInterface(k_RixMessages);

    initIO();

    if (!m_in.open(path)) {

		bbox[0] = 0;
		bbox[1] = 0;
		bbox[2] = 0;
		bbox[3] = 0;
		bbox[4] = 0;
		bbox[5] = 0;
        return;
    }
    // store the current to world xform matrix
    RtMatrix PrMt;
    RxTransform("current", "world", 1, PrMt);
    Mp[0][0] = PrMt[0][0]; Mp[0][1] = PrMt[0][1]; Mp[0][2] = PrMt[0][2]; Mp[0][3] = PrMt[0][3];
    Mp[1][0] = PrMt[1][0]; Mp[1][1] = PrMt[1][1]; Mp[1][2] = PrMt[1][2]; Mp[1][3] = PrMt[1][3];
    Mp[2][0] = PrMt[2][0]; Mp[2][1] = PrMt[2][1]; Mp[2][2] = PrMt[2][2]; Mp[2][3] = PrMt[2][3];
    Mp[3][0] = PrMt[3][0]; Mp[3][1] = PrMt[3][1]; Mp[3][2] = PrMt[3][2]; Mp[3][3] = PrMt[3][3];

    // loop over layers and partitions to access density
    m_in.getPartitionNames(m_partitions);
    const string dstr = "density";
    m_fieldName.push_back(dstr);
    
    BOOST_FOREACH (const string &partition, m_partitions) {

        if (!matchString(partition, m_fieldName)) {
            continue;
        }

        m_in.getScalarLayerNames(m_scalarLayers, partition);

        BOOST_FOREACH (const string &scalarLayer, m_scalarLayers) {

            if (!matchString(scalarLayer, m_fieldName)) {
                continue;
            }  
            Field<float>::Vec fScalarFields;
            fScalarFields = m_in.readScalarLayers<float>(partition, scalarLayer);

            if (fScalarFields.size() > 0) {

                ReadLayersAndSetupFields<float>(fScalarFields); 
            
                BOOST_FOREACH (Field<float>::Ptr float_field, fScalarFields) {
                    m_fbuffer = float_field;
                }
            }

            Field<double>::Vec dScalarFields;
            dScalarFields = m_in.readScalarLayers<double>(partition, scalarLayer);

            if (dScalarFields.size() > 0) {

                ReadLayersAndSetupFields<double>(dScalarFields); 

                BOOST_FOREACH (Field<double>::Ptr double_field, dScalarFields) {
                    m_dbuffer = double_field;
                }
            }

            Field<half>::Vec hScalarFields;
            hScalarFields = m_in.readScalarLayers<half>(partition, scalarLayer);

            if (hScalarFields.size() > 0) {

                ReadLayersAndSetupFields<half>(hScalarFields); 

                BOOST_FOREACH (Field<half>::Ptr half_field, hScalarFields) {
                    m_hbuffer = half_field;
                }
            }

        } 
    } 
#ifdef DEBUG
	msgs->Info("bbox: %f %f %f %f %f %f\n", 
				bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
#endif 
    // store velocity buffers
    vector<string> velFieldNames;
    vector<string> velFieldAttribs;
    vector<string> scalarLayers;

    vector<Field<float>::Ptr> vel_float_buf;
    vector<Field<double>::Ptr> vel_double_buf;
    vector<Field<half>::Ptr> vel_half_buf;
    // XXX: Test for vector field
    bool isVec = false;

    if (m_dataType == "float") { 

        vel_float_buf = VelocitySamples<float>(velFieldNames, velFieldAttribs, scalarLayers, isVec);

        if (vel_float_buf.size() > 1) {

            m_velx_mapping = vel_float_buf[0]->mapping().get();
            m_vely_mapping = vel_float_buf[1]->mapping().get();
            m_velz_mapping = vel_float_buf[2]->mapping().get();

            m_velx_fbuffer = vel_float_buf[0];
            m_vely_fbuffer = vel_float_buf[1];
            m_velz_fbuffer = vel_float_buf[2];

            collectVelocityBounds<float>(m_velx_fbuffer, m_vely_fbuffer, m_velz_fbuffer);
        } else { // it's a vector
            m_velx_mapping = vel_float_buf[0]->mapping().get();
            m_velx_fbuffer = vel_float_buf[0];
            // XXX: call func with T= V3f?
        }

        V3f vlen(m_vel_xmax, m_vel_ymax, m_vel_zmax);
        m_vel_length = vlen.length();
    }
    else if (m_dataType == "double") {

        vel_double_buf = VelocitySamples<double>(velFieldNames, velFieldAttribs, scalarLayers, isVec);

        m_velx_mapping = vel_double_buf[0]->mapping().get();
        m_vely_mapping = vel_double_buf[1]->mapping().get();
        m_velz_mapping = vel_double_buf[2]->mapping().get();

        m_velx_dbuffer = vel_double_buf[0];
        m_vely_dbuffer = vel_double_buf[1];
        m_velz_dbuffer = vel_double_buf[2];

        collectVelocityBounds<double>(m_velx_dbuffer, m_vely_dbuffer, m_velz_dbuffer);

        V3f vlen(m_vel_xmax, m_vel_ymax, m_vel_zmax);
        m_vel_length = vlen.length();
    }
    else if (m_dataType == "half") {

        vel_half_buf = VelocitySamples<half>(velFieldNames, velFieldAttribs, scalarLayers, isVec);

        m_velx_mapping = vel_half_buf[0]->mapping().get();
        m_vely_mapping = vel_half_buf[1]->mapping().get();
        m_velz_mapping = vel_half_buf[2]->mapping().get();

        m_velx_hbuffer = vel_half_buf[0];
        m_vely_hbuffer = vel_half_buf[1];
        m_velz_hbuffer = vel_half_buf[2];

        collectVelocityBounds<half>(m_velx_hbuffer, m_vely_hbuffer, m_velz_hbuffer);

        V3f vlen(m_vel_xmax, m_vel_ymax, m_vel_zmax);
        m_vel_length = vlen.length();
    }

}

template <typename T>
void F3DImplicitField::ReadLayersAndSetupFields(typename Field<T>::Vec sFields)
{
    BOOST_FOREACH (typename Field<T>::Ptr field, sFields) {

        /* Get info from cache */
        m_extents = field->extents();
        const V3i minres = m_extents.min;
        const V3i maxres = m_extents.max;
        minres.getValue(m_xmin, m_ymin, m_zmin);
        maxres.getValue(m_xmax, m_ymax, m_zmax);
        
        // data window
        m_datawin = field->dataWindow();
        const V3i dmaxres = m_datawin.max;
        dmaxres.getValue(m_dxmax, m_dymax, m_dzmax);
    
        m_mapping = field->mapping().get();

        // derive object-space bounds from extents
        V3d big((V3d)maxres + V3d(1));
        V3d maxvec = xformPoint(big, m_mapping, 0, 0);

        V3d small((V3d)minres + V3d(1));
        V3d minvec = xformPoint(small, m_mapping, 0, 0);
        // XXX: add padding in a more elegant fashion
        bbox[0] = minvec.x - 2;
        bbox[1] = maxvec.x + 2;
        bbox[2] = minvec.y - 2;
        bbox[3] = maxvec.y + 2;
        bbox[4] = minvec.z - 2;
        bbox[5] = maxvec.z + 2;
    
        m_dataType = field->dataTypeString();
        
    }
}

RtFloat F3DImplicitField::Eval(const RtPoint p)
{
    float ret;
    V3d wsP(p[0], p[1], p[2]);
    V3d vsP;
    
    m_mapping->worldToVoxel(wsP, vsP);
    //XXX: add security checks so out of bounds voxels are not sampled

    if (m_dataType == "float") {
        ret = m_finterpolator.sample(*m_fbuffer, vsP);
    }
    else if (m_dataType == "double") {
        ret = (float)m_dinterpolator.sample(*m_dbuffer, vsP);
    }
    else if (m_dataType == "half") {
        ret = (float)m_hinterpolator.sample(*m_hbuffer, vsP);
    }

    return ret;
}

void F3DImplicitField::EvalMultiple(int neval, float *result, int resultstride, const RtPoint *p) 
{
    for (int i = 0; i < neval; ++i) {
        *result = Eval(*p++);
        result += resultstride;
    }
}

void F3DImplicitField::Range(RtInterval val, const RtPoint corners[8],RtVolumeHandle h)
{
    val[0] = 0; val[1] = 0;

    bool ok = true;
    for(int i=0;i<8;i++)
    {
        ok = (ok &&(corners[i][0] < bbox[0]));
    }
    if(ok) return;

    for(int i=0;i<8;i++)
    {
        ok = (ok &&(corners[i][0] > bbox[1]));
    }
    if(ok) return;

    for(int i=0;i<8;i++)
    {
        ok = (ok &&(corners[i][1] < bbox[2]));
    }
    if(ok) return;

    for(int i=0;i<8;i++)
    {
        ok = (ok &&(corners[i][1] > bbox[3]));
    }
    if(ok) return;

    for(int i=0;i<8;i++)
    {
        ok = (ok &&(corners[i][2] < bbox[4]));
    }
    if(ok) return;

    for(int i=0;i<8;i++)
    {
        ok = (ok &&(corners[i][2] > bbox[5]));
    }
    if(ok) return;

    // DEF
    val[0]=-1e30;
    val[1]=1e30;
};

void F3DImplicitField::Motion(RtPoint result, const RtPoint p)
{
    V3d wsP(p[0], p[1], p[2]);
    V3d vsP, vsPx, vsPy, vsPz;
    V3f v_samps;
    float density;

    // transform to object space
    vsPx = xformPoint(wsP, m_velx_mapping, 0, 1);
    vsPy = xformPoint(wsP, m_vely_mapping, 0, 1);
    vsPz = xformPoint(wsP, m_velz_mapping, 0, 1);

    // sample buffers
    // XXX: Add security check so that out of bounds voxels are not sampled
    if (m_dataType == "float") {
        v_samps.x = m_finterpolator.sample(*m_velx_fbuffer, vsPx);
        v_samps.y = m_finterpolator.sample(*m_vely_fbuffer, vsPy);
        v_samps.z = m_finterpolator.sample(*m_velz_fbuffer, vsPz);
    } 
    else if (m_dataType == "double") {
        v_samps.x = m_dinterpolator.sample(*m_velx_dbuffer, vsPx);
        v_samps.y = m_dinterpolator.sample(*m_vely_dbuffer, vsPy);
        v_samps.z = m_dinterpolator.sample(*m_velz_dbuffer, vsPz);
    }
    else if (m_dataType == "half") {
        v_samps.x = m_hinterpolator.sample(*m_velx_hbuffer, vsPx);
        v_samps.y = m_hinterpolator.sample(*m_vely_hbuffer, vsPy);
        v_samps.z = m_hinterpolator.sample(*m_velz_hbuffer, vsPz);
    }

    result[0] = v_samps.x * m_blur;
    result[1] = v_samps.y * m_blur;
    result[2] = v_samps.z * m_blur;
}

void F3DImplicitField::BoxMotion(RtBound result, const RtBound b)
{
    V3d maxvel(m_vel_xmax, m_vel_ymax, m_vel_zmax);
    V3d vsPx = xformPoint(maxvel, m_velx_mapping, 0, 1);
    V3d vsPy = xformPoint(maxvel, m_vely_mapping, 0, 1);
    V3d vsPz = xformPoint(maxvel, m_velz_mapping, 0, 1);

    result[0] = b[0] - (vsPx.x * m_bbox_mod);
    result[1] = b[1] + (vsPx.x * m_bbox_mod);
    result[2] = b[2] - (vsPx.y * m_bbox_mod);
    result[3] = b[3] + (vsPx.y * m_bbox_mod);
    result[4] = b[4] - (vsPx.z * m_bbox_mod);
    result[5] = b[5] + (vsPx.z * m_bbox_mod);

#ifdef DEBUG
    msgs->Info("Motion bounds(%f %f %f %f %f %f)", result[0], result[1], result[2], result[3], result[4], result[5]);
#endif 
}

void F3DImplicitField::MotionMultiple(int neval, RtPoint *result, const RtPoint *p) 
{
    for (int i = 0; i < neval; ++i) {
        Motion(*result++, *p++);
    }
}

ImplicitVertexValue *F3DImplicitField::CreateVertexValue(const RtToken pname, int nvalue)
{
    if (m_dataType == "float") {
        return DeployField<float>(pname, nvalue);
    } 
    else if (m_dataType == "double") {
        return DeployField<double>(pname, nvalue);
    }
    else if (m_dataType == "half") {
        return DeployField<half>(pname, nvalue);
    }
}

// Utility function to transform points to world from current, and from voxel to world
V3d F3DImplicitField::xformPoint(V3d &Pt, FieldMapping *mapping, int frustrum, int fromobject)
{
    M44d vTw;
    if (frustrum != 0) {
        Field3D::FrustumFieldMapping::Ptr fmapping =
            boost::dynamic_pointer_cast<FrustumFieldMapping>(mapping);
        vTw = fmapping->cameraToWorld();
    } else {
        MatrixFieldMapping::Ptr xmapping = 
            boost::dynamic_pointer_cast<MatrixFieldMapping>(mapping);
        vTw = xmapping->voxelToWorld();
    }

    M44d oTw = Mp.inverse();

    M44d big;
    big.setTranslation(Pt);
    M44d mres;
    if (fromobject != 0) {
        mres.multiply(big, oTw, mres);
    } else {
        mres.multiply(big, vTw, mres);
    }
    V3d ret = mres.translation();

    return ret;
}

// Returns a buffer containing the velocity fields (X,Y,Z)
template <typename T>
vector<typename Field<T>::Ptr> F3DImplicitField::VelocitySamples(vector<string> &velFieldN, vector<string> &velAttrib, 
                                                                 vector<string> &sLayers, bool isVec)
{
    // query velocity fields
    const string name = "vel";
    size_t found;
    const string attribute_x = "x";
    const string attribute_y = "y";
    const string attribute_z = "z";

    // hold return buffer(s)
    vector<typename Field<T>::Ptr> ret;
    
    // velocity xyz buffers
    // XXX: Check for vector field containing all 3 velocities
    typedef typename Field<T>::Ptr SFieldList;
    typedef typename Field<FIELD3D_VEC3_T<T> >::Vec VFieldList;

    SFieldList sFields = m_in.readScalarLayers<T>(name);
    VFieldList vFields = m_in.readVectorLayers<T>(name);

    if (sFields.size() > 0) {
        isVec = false;
        typename Field<T>::Ptr vx = getField<T>(m_in, m_partitions, velFieldN, velAttrib, sLayers, attribute_x, name, isVec);
        typename Field<T>::Ptr vy = getField<T>(m_in, m_partitions, velFieldN, velAttrib, sLayers, attribute_y, name, isVec);
        typename Field<T>::Ptr vz = getField<T>(m_in, m_partitions, velFieldN, velAttrib, sLayers, attribute_z, name, isVec);
        
        ret.push_back(vx); ret.push_back(vy); ret.push_back(vz);
    } 
    else if (vFields.size() > 0) {
        isVec = true;
        // loop through vector fields
        for (typename VFieldList::const_iterator i = vFields.begin();
                i != vFields.end(); ++i) {
            // find velocity
            found = (**i).name.find(name);
            if (found != string::npos) {
                typename Field<T>::Ptr vel = getField<T>(m_in, m_partitions, velFieldN, velAttrib, sLayers, 
                                                        (**i).attribute, name, isVec);
            }
        }
        ret.push_back(vel);
    }
    // might contain 1 or 3 items
    return ret;
}

template <typename T>
VertexField<T>* F3DImplicitField::DeployField(const RtToken &pname, int &nvalue)
{
    LinearFieldInterp<T> interp;
    const string fuel = "fuel";
    const string temp = "temperature";
    const string burn = "burn";
    const string heat = "heat";

    // only test nvalue for scalars
    if (nvalue == 1 && !strcmp(pname, "float fuel")) {
        return new VertexField<T>(m_in, m_mapping, interp, m_partitions, Mp, fuel);
    }
    else if (nvalue == 1 && !strcmp(pname, "float temperature")) {
        return new VertexField<T>(m_in, m_mapping, interp, m_partitions, Mp, temp);
    }
    else if (nvalue == 1 && !strcmp(pname, "float heat")) {
        return new VertexField<T>(m_in, m_mapping, interp, m_partitions, Mp, heat);
    }
    else if (nvalue == 1 && !strcmp(pname, "float burn")) {
        return new VertexField<T>(m_in, m_mapping, interp, m_partitions, Mp, burn);
    } else {
        return NULL;
    }
}

bool matchString(const string &str, const vector<string> &patterns)
{
  // If patterns is empty all strings match
  if (patterns.size() == 0) {
    return true;
  }
  // Check all patterns
  BOOST_FOREACH (const string &pattern, patterns) {
    if (fnmatch(pattern.c_str(), str.c_str(), 0) != FNM_NOMATCH) {
      return true;
    }
  }
  // If no pattern matched return false
  return false;
}

// Non member function to query a named scalar field
template <typename T>
typename Field<T>::Ptr getField(Field3DInputFile &file, vector<string> &partitionsVec, vector<string> &fieldName,
        vector<string> &attribName, vector<string> &scalarLayers, const string &attrib, const string &name, bool isVec)
{
    file.getPartitionNames(partitionsVec);

    attribName.push_back(attrib);
    fieldName.push_back(name);

    BOOST_FOREACH (const string &partition, partitionsVec) {

        if (!matchString(partition, fieldName)) {
            continue;
        }

        file.getScalarLayerNames(scalarLayers, partition);

        BOOST_FOREACH (const string &scalarLayer, scalarLayers) {

            if (!matchString(scalarLayer, attribName)) {
                continue;
            }

            if (isVec) {
                typename Field<FIELD3D_VEC3_T<T> >::Vec vectorFields;
                vectorFields = file.readVectorLayers<T>(partition, scalarLayer);

                if (vectorFields.size() > 0) {

                    BOOST_FOREACH (typename Field<T>::Ptr vfield, vectorFields) {

                        return vfield;
                    }
                }
            } else {
                typename Field<T>::Vec fScalarFields;
                fScalarFields = file.readScalarLayers<T>(partition, scalarLayer);

                if (fScalarFields.size() > 0) {

                    BOOST_FOREACH (typename Field<T>::Ptr sfield, fScalarFields) {

                        return sfield;
                    }
                }
            }

        } // layers
    } // partitions
    return NULL;
}

// XXX: T should be V3f, V3d, or V3h(?)
template <typename T>
void F3DImplicitField::collectVectorVelocityBounds(typename Field<T>::Ptr vel)
{
    size_t iX, iY, iZ;
    T vx;
    vector<T> velocities;

    // iterate voxels
    // XXX: Using the bounds from density may cause out of bounds calcs
    for(iZ = 0; iZ < m_zmax; iZ++)
    {
      for (iY = 0; iY < m_ymax; iY++)
      {
          for (iX = 0; iX < m_xmax; iX++)
          {
              // sample the values
              vx = vel->value(iX, iY, iZ);
              // push the velocities into vectors (one component from each V3d)
              velocities.push_back(vx);
          }
      }
    }

    vector<T>::const_iterator i = velocities.begin();
    for (; i != velocities.end(); ++i) {
        //XXX: figure out max value of a V3f
        //m_vel_xmax = *max_element(*
    }
    // find the max value in each axis, and store it in member data
    //m_vel_xmax = *max_element(velocitiesx.begin(), velocitiesx.end());
    //m_vel_ymax = *max_element(velocitiesy.begin(), velocitiesy.end());
    //m_vel_zmax = *max_element(velocitiesz.begin(), velocitiesz.end());
}

template <typename T>
void F3DImplicitField::collectVelocityBounds(typename Field<T>::Ptr xbuf,
                                             typename Field<T>::Ptr ybuf,
                                             typename Field<T>::Ptr zbuf)
{
    size_t iX, iY, iZ;
    T vx, vy, vz;
    vector<T> velocitiesx;
    vector<T> velocitiesy;
    vector<T> velocitiesz;

    // iterate voxels
    // XXX: Using the bounds from density may cause out of bounds calcs
    for(iZ = 0; iZ < m_zmax; iZ++)
    {
      for (iY = 0; iY < m_ymax; iY++)
      {
          for (iX = 0; iX < m_xmax; iX++)
          {
              // sample the values
              vx = xbuf->value(iX, iY, iZ);
              vy = ybuf->value(iX, iY, iZ);
              vz = zbuf->value(iX, iY, iZ);
              // push the velocities into vectors (one component from each V3d)
              velocitiesx.push_back(vx);
              velocitiesy.push_back(vy);
              velocitiesz.push_back(vz);
          }
      }
    }
    // find the max value in each axis, and store it in member data
    m_vel_xmax = *max_element(velocitiesx.begin(), velocitiesx.end());
    m_vel_ymax = *max_element(velocitiesy.begin(), velocitiesy.end());
    m_vel_zmax = *max_element(velocitiesz.begin(), velocitiesz.end());
}

FIELDCREATE
{
    float blur = 1;
    float bbox_mod = 1;
    if(nfloat > 0) blur = float0[0];
    if(nfloat > 1) bbox_mod = float0[1];
	if(nstring > 0)
		return new F3DImplicitField(string[0], blur, bbox_mod);
	else
		return 0;
}
