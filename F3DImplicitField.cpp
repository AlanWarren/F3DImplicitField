#include <iostream>
#include <algorithm>
#include <cstdio>
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

#include <OpenEXR/ImathVec.h>
#include <OpenEXR/ImathBox.h>
#include <OpenEXR/ImathBoxAlgo.h>

#include "ImplicitField.h"
#include "rx.h"
#include "RixInterfaces.h"

using namespace std;
using namespace Field3D;

bool matchString(const string&, const vector<string>&);

template <typename T>
typename Field<T>::Ptr getField(Field3DInputFile&, vector<string>&, vector<string>&, vector<string>&, 
        vector<string>&, const string&, const string&);

template <typename T>
typename Field<FIELD3D_VEC3_T<T> >::Ptr getVectorField(Field3DInputFile&, vector<string>&, vector<string>&, vector<string>&, 
        vector<string>&, const string&, const string&);

/* VertexField class */
template <class T>
class VertexField : public ImplicitVertexValue
{
    public:
        VertexField(Field3DInputFile&, FieldMapping*, LinearFieldInterp<T>, vector<string>&, M44d, const string&);
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

// constructor initializes data necessary to query any named field (fuel, heat, temperature, etc.)
template <class T>
VertexField<T>::VertexField(Field3DInputFile &file, FieldMapping *mapping, LinearFieldInterp<T> interp,
        vector<string> &partitionsVec, M44d oTw, const string &fieldName) : 
    t_mapping(mapping), 
    t_interpolator(interp), 
    t_partitions(partitionsVec),
    t_oTw(oTw),
    t_field(fieldName)
{
    t_buffer = getField<T>(file, t_partitions, t_fieldName, t_attribName, t_scalarLayers, t_field, t_field);
}

template <class T>
RtFloat VertexField<T>::EvalField(const RtPoint p) 
{
    float ret;
    V3d wsP(p[0], p[1], p[2]);
    V3d vsP;

    t_mapping->worldToVoxel(wsP, vsP);
    return t_interpolator.sample(*t_buffer, vsP);
}

class F3DImplicitField : public ImplicitField
{
    public:
    	F3DImplicitField(string path, string field, float blur, float bbox_mod, float blur_cubic, float field_cubic);
    	~F3DImplicitField(){}
    
    	virtual RtFloat Eval(const RtPoint p);
        virtual void EvalMultiple(int neval, float *result, int resultstride, const RtPoint *p);
        virtual void GradientEval(RtPoint result, const RtPoint p) {result[0] = 0; result[1] = 0; result[2] = 0;}
    	virtual void Range(RtInterval result, const RtPoint corners[8], const RtVolumeHandle h);
        virtual void Motion(RtPoint result, const RtPoint p);
        virtual void BoxMotion(RtBound result, const RtBound b);
        virtual void MotionMultiple(int neval, RtPoint *result, const RtPoint *p);
        ImplicitVertexValue *CreateVertexValue(const RtToken name, int nvalue);
   
        V3d xformPoint(V3d &Pt, FieldMapping *mapping, int fromobject, int fromworld);

        template <typename T>
        void ReadLayersAndSetupFields(typename Field<T>::Vec sFields);

        template <typename T>
        typename Field<FIELD3D_VEC3_T<T> >::Ptr SetupVectorVelocityFields(vector<string> &velFieldNames,
                                                                          vector<string> &velFieldAttribs,
                                                                          vector<string> &scalarLayers,
                                                                          const string   &velpartition);

        template <typename T>
        vector<typename Field<T>::Ptr> SetupScalarVelocityFields(vector<string> &velFieldN, 
                                                                 vector<string> &velAttrib, 
                                                                 vector<string> &partVec,
                                                                 vector<string> &sLayers);

        template <typename T>
        VertexField<T>* NamedVertexField(const RtToken &pname, int &nvalue);

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
        // argv
        float m_blur; 
        float m_bbox_mod;
        float m_blur_cubic;
        float m_field_cubic;
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
        CubicFieldInterp<float> m_cubicfinterpolator;
        CubicFieldInterp<double> m_cubicdinterpolator;
        CubicFieldInterp<half> m_cubichinterpolator;
        // density buffers
        Field<float>::Ptr m_fbuffer;
        Field<double>::Ptr m_dbuffer;
        Field<half>::Ptr m_hbuffer;
        // vertex fields
        vector<string> m_vertexfields;
        // scalar velocity buffers
        Field<float>::Ptr m_velx_fbuffer;
        Field<float>::Ptr m_vely_fbuffer;
        Field<float>::Ptr m_velz_fbuffer;
        Field<double>::Ptr m_velx_dbuffer;
        Field<double>::Ptr m_vely_dbuffer;
        Field<double>::Ptr m_velz_dbuffer;
        Field<half>::Ptr m_velx_hbuffer;
        Field<half>::Ptr m_vely_hbuffer;
        Field<half>::Ptr m_velz_hbuffer;
        // vector velocity buffers
        Field<V3f>::Ptr m_vel_fbuffer;
        Field<V3d>::Ptr m_vel_dbuffer;
        Field<V3h>::Ptr m_vel_hbuffer;
        bool m_vector_field;
        // velocity mappings
        FieldMapping *m_vel_mapping;
        // scalar mappings
        FieldMapping *m_velx_mapping;
        FieldMapping *m_vely_mapping;
        FieldMapping *m_velz_mapping;
        // bounds modifiers
        float m_vel_xmax, m_vel_ymax, m_vel_zmax;
        float m_vel_xmin, m_vel_ymin, m_vel_zmin;
        double m_vel_length;
};

F3DImplicitField::F3DImplicitField(string path, string field, float blur, float bbox_mod, float blur_cubic, float field_cubic) : 
    m_blur(blur), m_bbox_mod(bbox_mod), m_blur_cubic(blur_cubic), m_field_cubic(field_cubic)
{
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

    if (!field.empty()) {
        m_fieldName.push_back(field);
    } else {
        m_fieldName.push_back(dstr);
    }

    // loop over once to gather all field names
    BOOST_FOREACH (const string &part, m_partitions) {
        m_vertexfields.push_back(part);
    }

    // Density Field
    BOOST_FOREACH (const string &partition, m_partitions) {

        if (!matchString(partition, m_fieldName)) {
            continue;
        }

        m_in.getScalarLayerNames(m_scalarLayers, partition);

        BOOST_FOREACH (const string &scalarLayer, m_scalarLayers) {

            if (!matchString(scalarLayer, m_fieldName)) {
                continue;
            }  
            // if our f3d file contains float data
            Field<float>::Vec fScalarFields;
            fScalarFields = m_in.readScalarLayers<float>(partition, scalarLayer);

            if (fScalarFields.size() > 0) {

                ReadLayersAndSetupFields<float>(fScalarFields); 
            
                BOOST_FOREACH (Field<float>::Ptr float_field, fScalarFields) {
                    m_fbuffer = float_field;
                }
            }

            // if our f3d file contains double data
            Field<double>::Vec dScalarFields;
            dScalarFields = m_in.readScalarLayers<double>(partition, scalarLayer);

            if (dScalarFields.size() > 0) {

                ReadLayersAndSetupFields<double>(dScalarFields); 

                BOOST_FOREACH (Field<double>::Ptr double_field, dScalarFields) {
                    m_dbuffer = double_field;
                }
            }

            // if our f3d file contains double data
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
   
    // Velocity Field(s)
    m_vector_field = false;
    vector<string> velFieldNames;
    vector<string> velFieldAttribs;
    vector<string> scalarLayers;
    vector<string> partVec;
    const string   velpartition = "vel";
    const string   velattrib = "x";

    if (m_dataType == "float") { 

        Field<V3f>::Vec vFields;
        vFields = m_in.readVectorLayers<float>(velpartition, velpartition);

        Field<float>::Vec sFields;
        sFields = m_in.readScalarLayers<float>(velpartition, velattrib);

        if (vFields.size() > 0) {
            m_vector_field = true;
            m_vel_fbuffer = SetupVectorVelocityFields<float>(velFieldNames,
                                                             velFieldAttribs,
                                                             scalarLayers,
                                                             velpartition);
        } else if (sFields.size() > 0) {
            vector<Field<float>::Ptr> vecbuf;
            vecbuf = SetupScalarVelocityFields<float>(velFieldNames, 
                                                      velFieldAttribs, 
                                                      partVec,
                                                      scalarLayers);
            m_velx_fbuffer = vecbuf[0];
            m_vely_fbuffer = vecbuf[1];
            m_velz_fbuffer = vecbuf[2];
        }
    } else if (m_dataType == "double") {
        
        Field<V3d>::Vec vFields;
        vFields = m_in.readVectorLayers<double>(velpartition, velpartition);

        Field<double>::Vec sFields;
        sFields = m_in.readScalarLayers<double>(velpartition, velattrib);

        if (vFields.size() > 0) {
            m_vector_field = true;
            m_vel_dbuffer = SetupVectorVelocityFields<double>(velFieldNames,
                                                              velFieldAttribs,
                                                              scalarLayers,
                                                              velpartition);
        } else if (sFields.size() > 0) {
            vector<Field<double>::Ptr> vecbuf;
            vecbuf = SetupScalarVelocityFields<double>(velFieldNames, 
                                                       velFieldAttribs, 
                                                       partVec,
                                                       scalarLayers);
            m_velx_dbuffer = vecbuf[0];
            m_vely_dbuffer = vecbuf[1];
            m_velz_dbuffer = vecbuf[2];
        }
    } else if (m_dataType == "half") {

        Field<V3h>::Vec vFields;
        vFields = m_in.readVectorLayers<half>(velpartition, velpartition);

        Field<half>::Vec sFields;
        sFields = m_in.readScalarLayers<half>(velpartition, velattrib);

        if (vFields.size() > 0) {
            m_vector_field = true;
            m_vel_hbuffer = SetupVectorVelocityFields<half>(velFieldNames,
                                                            velFieldAttribs,
                                                            scalarLayers,
                                                            velpartition);
        } else if (sFields.size() > 0) {
            vector<Field<half>::Ptr> vecbuf;
            vecbuf = SetupScalarVelocityFields<half>(velFieldNames, 
                                                     velFieldAttribs, 
                                                     partVec,
                                                     scalarLayers);
            m_velx_hbuffer = vecbuf[0];
            m_vely_hbuffer = vecbuf[1];
            m_velz_hbuffer = vecbuf[2];
        }
    }
    // bounds
	msgs->Info("Bounding Box: %f %f %f %f %f %f\n", 
				bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
}

void F3DImplicitField::Motion(RtPoint result, const RtPoint p)
{
    V3d wsP(p[0], p[1], p[2]);
    //RtPoint mp[1] = {.0f, .0f, .0f};
    //RxTransformPoints("object", "world", 1 p, 0.0f);
    V3d vsP, vsPx, vsPy, vsPz;
    V3f ret;

    const string   velpartition = "vel";
    const string   velattrib = "x";

    // Define mappings 
    if (m_vector_field) {

        m_vel_mapping->worldToVoxel(wsP, vsP);

    } else {
        m_velx_mapping->worldToVoxel(wsP, vsPx);
        m_vely_mapping->worldToVoxel(wsP, vsPy);
        m_velz_mapping->worldToVoxel(wsP, vsPz);
    }

    // Sample velocity field(s)
    if (m_dataType == "float") {
        LinearFieldInterp<V3f> f3_linearinterpolator;
        CubicFieldInterp<V3f> f3_cubicinterpolator;

        V3f v_samps;

        if (m_vector_field) {
            // cubic or linear interpolation
            if (m_blur_cubic != 0) {
                v_samps = f3_cubicinterpolator.sample(*m_vel_fbuffer, vsP);
            } else {
                v_samps = f3_linearinterpolator.sample(*m_vel_fbuffer, vsP);
            }

        } else {

            if (m_blur_cubic != 0) {
                v_samps.x = m_cubicfinterpolator.sample(*m_velx_fbuffer, vsPx);
                v_samps.y = m_cubicfinterpolator.sample(*m_vely_fbuffer, vsPy);
                v_samps.z = m_cubicfinterpolator.sample(*m_velz_fbuffer, vsPz);
            } else {
                v_samps.x = m_finterpolator.sample(*m_velx_fbuffer, vsPx);
                v_samps.y = m_finterpolator.sample(*m_vely_fbuffer, vsPy);
                v_samps.z = m_finterpolator.sample(*m_velz_fbuffer, vsPz);
            }
        }

        ret.x = v_samps.x;
        ret.y = v_samps.y;
        ret.z = v_samps.z;

    } else if (m_dataType == "double") {
        LinearFieldInterp<V3d> d3_linearinterpolator;
        CubicFieldInterp<V3d> d3_cubicinterpolator;

        V3d v_samps;

        if (m_vector_field) {

            if (m_blur_cubic != 0) {
                v_samps = d3_cubicinterpolator.sample(*m_vel_dbuffer, vsP);
            } else {
                v_samps = d3_linearinterpolator.sample(*m_vel_dbuffer, vsP);
            }

        } else {

            if (m_blur_cubic != 0) {
                v_samps.x = m_cubicdinterpolator.sample(*m_velx_dbuffer, vsPx);
                v_samps.y = m_cubicdinterpolator.sample(*m_vely_dbuffer, vsPy);
                v_samps.z = m_cubicdinterpolator.sample(*m_velz_dbuffer, vsPz);
            } else {
                v_samps.x = m_dinterpolator.sample(*m_velx_dbuffer, vsPx);
                v_samps.y = m_dinterpolator.sample(*m_vely_dbuffer, vsPy);
                v_samps.z = m_dinterpolator.sample(*m_velz_dbuffer, vsPz);
            }
        }

        ret.x = static_cast<float>(v_samps.x);
        ret.y = static_cast<float>(v_samps.y);
        ret.z = static_cast<float>(v_samps.z);

    } else if (m_dataType == "half") {
        LinearFieldInterp<V3h> h3_linearinterpolator;
        CubicFieldInterp<V3h> h3_cubicinterpolator;

        V3h v_samps;

        if (m_vector_field) {

            if (m_blur_cubic != 0) {
                v_samps = h3_cubicinterpolator.sample(*m_vel_hbuffer, vsP);
            } else {
                v_samps = h3_linearinterpolator.sample(*m_vel_hbuffer, vsP);
            }

        } else {

            if (m_blur_cubic != 0) {
                v_samps.x = m_cubichinterpolator.sample(*m_velx_hbuffer, vsPx);
                v_samps.y = m_cubichinterpolator.sample(*m_vely_hbuffer, vsPy);
                v_samps.z = m_cubichinterpolator.sample(*m_velz_hbuffer, vsPz);
            } else {
                v_samps.x = m_hinterpolator.sample(*m_velx_hbuffer, vsPx);
                v_samps.y = m_hinterpolator.sample(*m_vely_hbuffer, vsPy);
                v_samps.z = m_hinterpolator.sample(*m_velz_hbuffer, vsPz);
            }
        }

        ret.x = static_cast<float>(v_samps.x);
        ret.y = static_cast<float>(v_samps.y);
        ret.z = static_cast<float>(v_samps.z);
    }

	ret.x /= 24.0/0.3;
	ret.y /= 24.0/0.3;
	ret.z /= 24.0/0.3;

    result[0] = ret.x * m_blur;
    result[1] = ret.y * m_blur;
    result[2] = ret.z * m_blur;

}

void F3DImplicitField::BoxMotion(RtBound result, const RtBound b)
{
    // Take the difference between our original bbox and motion blurred bbox
    V3f minlen(m_vel_xmin, m_vel_ymin, m_vel_zmin);
    V3f maxlen(m_vel_xmax, m_vel_ymax, m_vel_zmax);
    float minlenf = minlen.length();
    float maxlenf = maxlen.length();

    V3d diffmax(m_vel_xmax - bbox[1], m_vel_ymax - bbox[3], m_vel_zmax - bbox[5]);
    V3d diffmin(m_vel_xmin - bbox[0], m_vel_ymin - bbox[2], m_vel_zmin - bbox[4]);

    // Offset by length seems to work the best..
    result[0] = b[0] - (maxlenf * m_bbox_mod);
    result[1] = b[1] + (maxlenf * m_bbox_mod);
    result[2] = b[2] - (maxlenf * m_bbox_mod);
    result[3] = b[3] + (maxlenf * m_bbox_mod);
    result[4] = b[4] - (maxlenf * m_bbox_mod);
    result[5] = b[5] + (maxlenf * m_bbox_mod);
    
    //msgs->Info("Motion Bounds (%f %f %f %f %f %f)", result[0], result[1], result[2], result[3], result[4], result[5]);
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
        return NamedVertexField<float>(pname, nvalue);
    } 
    else if (m_dataType == "double") {
        return NamedVertexField<double>(pname, nvalue);
    }
    else if (m_dataType == "half") {
        return NamedVertexField<half>(pname, nvalue);
    }
}

// Returns fields both defined in RIB, and present in the cache (besides vel & density)
template <typename T>
VertexField<T>* F3DImplicitField::NamedVertexField(const RtToken &pname, int &nvalue)
{
    string flt = "float ";
    string cmp, nm;
    LinearFieldInterp<T> interp;
    // iterate through all fields present in cache
    vector<string>::const_iterator vit = m_vertexfields.begin();
    for (; vit != m_vertexfields.end(); ++vit)
    {
        if (*vit != "vel" && *vit != "density")
        { 
            // &pname contains "float fieldname", so I compose a string to compare
            cmp = flt + *vit;
            nm = *vit;
            // If the field is scalar AND present in our RIB declaration
            if (nvalue == 1 && !strcmp(pname, cmp.c_str())) {
                return new VertexField<T>(m_in, m_mapping, interp, m_partitions, Mp, nm);
            } 
        }
    }
    return NULL;
}

// Setup extents and bounds
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

        // convert bounds to world space
        V3d big((V3d)maxres + V3d(1));
        V3d maxvec = xformPoint(big, m_mapping, 0, 0);

        V3d small((V3d)minres + V3d(1));
        V3d minvec = xformPoint(small, m_mapping, 0, 0);

        bbox[0] = minvec.x;
        bbox[1] = maxvec.x;
        bbox[2] = minvec.y;
        bbox[3] = maxvec.y;
        bbox[4] = minvec.z;
        bbox[5] = maxvec.z;
    
        m_dataType = field->dataTypeString();
        
    }
}

RtFloat F3DImplicitField::Eval(const RtPoint p)
{
    float ret;
    V3d wsP(p[0], p[1], p[2]);
    V3d vsP;

    m_mapping->worldToVoxel(wsP, vsP);

    if (m_dataType == "float") {

        if (m_field_cubic != 0) {
            ret = m_cubicfinterpolator.sample(*m_fbuffer, vsP);
        } else {
            ret = m_finterpolator.sample(*m_fbuffer, vsP);
        }

    }
    else if (m_dataType == "double") {

        if (m_field_cubic != 0) {
            ret = (float)m_cubicdinterpolator.sample(*m_dbuffer, vsP);
        } else {
            ret = (float)m_dinterpolator.sample(*m_dbuffer, vsP);
        }

    }
    else if (m_dataType == "half") {

        if (m_field_cubic != 0) {
            ret = (float)m_cubichinterpolator.sample(*m_hbuffer, vsP);
        } else {
            ret = (float)m_hinterpolator.sample(*m_hbuffer, vsP);
        }

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

/* returns the vector velocity buffer and sets member
 * data for min / max velocities sampled from buffer,
 * as well as the max length of said values */
template <typename T>
typename Field<FIELD3D_VEC3_T<T> >::Ptr 
F3DImplicitField::SetupVectorVelocityFields(vector<string> &velFieldNames,
                                            vector<string> &velFieldAttribs,
                                            vector<string> &scalarLayers,
                                            const string   &velpartition)
{
    typename Field<FIELD3D_VEC3_T<T> >::Ptr buffer;
    buffer = getVectorField<T>(m_in, m_partitions, velFieldNames,
                               velFieldAttribs, scalarLayers, 
                               velpartition, velpartition);

    m_vel_mapping = buffer->mapping().get();
    Box3i xext = buffer->extents();
    const V3i xvx = xext.max;

    size_t iX, iY, iZ;
    vector<T> velx, vely, velz;
    FIELD3D_VEC3_T<T> bufVel;
    vector<FIELD3D_VEC3_T<T> > velocities;

    // iterate voxels and store velocity data from cache
    for(iZ = 0; iZ < xvx.z; iZ++)
    {
        for (iY = 0; iY < xvx.y; iY++)
        {
            for (iX = 0; iX < xvx.x; iX++)
            {
                // sample the values
                bufVel = buffer->value(iX, iY, iZ);
                velocities.push_back(bufVel);
            }
        }
    }

    // loop over V3f and store x y z in separate vectors
    typename vector<FIELD3D_VEC3_T<T> >::const_iterator vit = velocities.begin();
    for (; vit != velocities.end(); ++vit) 
    {
        velx.push_back(vit->x);
        vely.push_back(vit->y);
        velz.push_back(vit->z);
    }

    // find the max value
    m_vel_xmax = *max_element(velx.begin(), velx.end());
    m_vel_ymax = *max_element(vely.begin(), vely.end());
    m_vel_zmax = *max_element(velz.begin(), velz.end());
    // find the min value
    m_vel_xmin = *min_element(velx.begin(), velx.end());
    m_vel_ymin = *min_element(vely.begin(), vely.end());
    m_vel_zmin = *min_element(velz.begin(), velz.end());

    // get the length of our max velocity vec
    V3d vlen(m_vel_xmax, m_vel_ymax, m_vel_zmax);
    m_vel_length = vlen.length();

    return buffer;

}

template <typename T>
vector<typename Field<T>::Ptr> 
F3DImplicitField::SetupScalarVelocityFields(vector<string> &velFieldN, 
                                            vector<string> &velAttrib, 
                                            vector<string> &partVec,
                                            vector<string> &sLayers)
{
    vector<typename Field<T>::Ptr> vecbuf;

    const string name = "vel";
    const string attribute_x = "x";
    const string attribute_y = "y";
    const string attribute_z = "z";
    
    // velocity xyz buffers
    typename Field<T>::Ptr vx = getField<T>(m_in, partVec, velFieldN, velAttrib, sLayers, attribute_x, name);
    typename Field<T>::Ptr vy = getField<T>(m_in, partVec, velFieldN, velAttrib, sLayers, attribute_y, name);
    typename Field<T>::Ptr vz = getField<T>(m_in, partVec, velFieldN, velAttrib, sLayers, attribute_z, name);
    // store buffers for return
    vecbuf.push_back(vx); vecbuf.push_back(vy); vecbuf.push_back(vz);

    m_velx_mapping = vx->mapping().get();
    m_vely_mapping = vy->mapping().get();
    m_velz_mapping = vz->mapping().get();

    size_t iX, iY, iZ;
    size_t jX, jY, jZ;
    size_t kX, kY, kZ;
    T velx, vely, velz;
    vector<T> velocitiesx;
    vector<T> velocitiesy;
    vector<T> velocitiesz;
    Box3i xExt, yExt, zExt;
    xExt = vx->extents(); yExt = vy->extents(); zExt = vz->extents();
    const V3i xMaxres = xExt.max;
    const V3i yMaxres = yExt.max;
    const V3i zMaxres = zExt.max;

    for(iZ = 0; iZ < xMaxres.z; iZ++)
    {
        for (iY = 0; iY < xMaxres.y; iY++)
        {
            for (iX = 0; iX < xMaxres.x; iX++)
            {
                // sample the values
                velx = vx->value(iX, iY, iZ);
                vely = vy->value(iX, iY, iZ);
                velz = vz->value(iX, iY, iZ);

                velocitiesx.push_back(velx);
                velocitiesy.push_back(vely);
                velocitiesz.push_back(velz);
            }
        }
    }

    // find the max value in each axis, and store it in member data
    m_vel_xmax = *max_element(velocitiesx.begin(), velocitiesx.end());
    m_vel_ymax = *max_element(velocitiesy.begin(), velocitiesy.end());
    m_vel_zmax = *max_element(velocitiesz.begin(), velocitiesz.end());
    // find the min value
    m_vel_xmin = *min_element(velocitiesx.begin(), velocitiesx.end());
    m_vel_ymin = *min_element(velocitiesy.begin(), velocitiesy.end());
    m_vel_zmin = *min_element(velocitiesz.begin(), velocitiesz.end());

    // acquire the length and store it
    V3f vlen(m_vel_xmax, m_vel_ymax, m_vel_zmax);
    m_vel_length = vlen.length();

    return vecbuf;
}

template <typename T>
typename Field<T>::Ptr getField(Field3DInputFile &file, vector<string> &partitionsVec, vector<string> &fieldName,
        vector<string> &attribName, vector<string> &scalarLayers, const string &attrib, const string &name)
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

            typename Field<T>::Vec fScalarFields;
            fScalarFields = file.readScalarLayers<T>(partition, scalarLayer);

            if (fScalarFields.size() > 0) {

                BOOST_FOREACH (typename Field<T>::Ptr field, fScalarFields) {

                    return field;
                }
            }

        } 
    } 
    return NULL;
}

/* returns a pointer to a velocity buffer (XYZ) */
template <typename T>
typename Field<FIELD3D_VEC3_T<T> >::Ptr getVectorField(Field3DInputFile &file, vector<string> &partitionsVec, vector<string> &fieldName,
        vector<string> &attribName, vector<string> &vectorLayers, const string &attrib, const string &name)
{
    file.getPartitionNames(partitionsVec);

    attribName.push_back(attrib);
    fieldName.push_back(name);

    BOOST_FOREACH (const string &partition, partitionsVec) {

        if (!matchString(partition, fieldName)) {
            continue;
        }

        file.getVectorLayerNames(vectorLayers, partition);

        BOOST_FOREACH (const string &vectorLayer, vectorLayers) {

            if (!matchString(vectorLayer, attribName)) {
                continue;
            }  

            typename Field<FIELD3D_VEC3_T<T> >::Vec fVectorFields;
            fVectorFields = file.readVectorLayers<T>(partition, vectorLayer);

            if (fVectorFields.size() > 0) {

                BOOST_FOREACH (typename Field<FIELD3D_VEC3_T<T> >::Ptr field, fVectorFields) {

                    return field;
                }
            }

        } 
    } 
    return NULL;
}


V3d F3DImplicitField::xformPoint(V3d &Pt, FieldMapping *mapping, int fromobject, int fromworld)
{

    MatrixFieldMapping::Ptr xmapping = 
        boost::dynamic_pointer_cast<MatrixFieldMapping>(mapping);

    M44d vTw = xmapping->voxelToWorld();
    M44d wTo = Mp.inverse();

    M44d big;
    M44d mres;
    // matrix containing our point
    big.makeIdentity();
    big.setTranslation(Pt);

    if (fromobject != 0) {
        // world to object transform
        mres.multiply(big, wTo, mres);
    } else if (fromworld != 0) {
        // object to world transform
        mres.multiply(big, Mp, mres);
    } else {
        // voxel to world transform
        mres.multiply(big, vTw, mres);
    }
    V3d ret = mres.translation();

    return ret;
}

FIELDCREATE
{
    float blur = 1;
    float bbox_mod = 1;
    float blur_cubic = 0;
    float field_cubic = 0;
    if(nfloat > 0) blur = float0[0];
    if(nfloat > 1) bbox_mod = float0[1];
    if(nfloat > 2) blur_cubic = float0[2];
    if(nfloat > 3) field_cubic = float0[3];
	if(nstring > 0)
		return new F3DImplicitField(string[0], string[1], blur, bbox_mod, blur_cubic, field_cubic);
	else
		return 0;
}

