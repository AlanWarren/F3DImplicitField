/* f3d2prman
 * author: Alan Warren <email@alan-warren.com>
 * comments: This is a runprogram for F3DImplicitField
 */

#include <iostream>
#include <string>

#include <fnmatch.h>

#include <boost/pointer_cast.hpp>
#include <boost/foreach.hpp>

#include <Field3D/DenseField.h>
#include <Field3D/MACField.h>
#include <Field3D/SparseField.h>
#include <Field3D/InitIO.h>
#include <Field3D/Field3DFile.h>

#include "rx.h"
#include "RixInterfaces.h"

using namespace std;
using namespace Field3D;

bool matchString(const string&, const vector<string>&);

template <typename T>
void ReadLayersAndSetupBounds(typename Field<T>::Vec, V3d&, V3d&, FieldMapping*);

V3d xformPoint(V3d&, FieldMapping*);

int main(int argc, char **argv)
{
    char buf[4096];
    float detail;
    int Mode;

    Field3D::initIO();

	RixContext *rixCtx = RxGetRixContext();
	RixMessages *msgs = (RixMessages*)rixCtx->GetRixInterface(k_RixMessages);

    string filename;
    float blur, bbox_mod, threshold, blur_cubic, field_cubic;
    int doblur = 0;
    float shutter_open = 0;
    float shutter_close = 1;
    int blur_distance = 2;

    if (argc == 11) {
      filename = string(argv[1]);
      blur = atof(argv[2]);
      bbox_mod = atof(argv[3]);
      blur_cubic = atof(argv[4]);
      field_cubic = atof(argv[5]);
      threshold = atof(argv[6]);
      doblur  = atoi(argv[7]);
      shutter_open = atof(argv[8]);
      shutter_close = atof(argv[9]);
      blur_distance = atoi(argv[10]);
    } else {
        msgs->Error("Not enough arguments specified");
        msgs->Info("Usage: f3d2prman filename blur bbox_mod blur_cubic field_cubic threshold doblur shutter_open shutter_close blur_distance");
        return 1;
    }

    Field3DInputFile in;

    if (!in.open(filename)) {
        msgs->Error("Couldn't open %s", filename.c_str());
        return 1;
    }
    // Prepare for field processing ---
    const string fieldtype = "density";
    vector<string> fieldName;
    vector<string> attribName;
    vector<string> layers;
    vector<string> partitions;
    in.getPartitionNames(partitions);
    FieldMapping *mapping;
    V3d minBound;
    V3d maxBound;

    // vertex fields
    vector<string> vertexFields;

    while(!std::cin.eof()) {

        sscanf(buf, "%f %d", &detail, &Mode);
        // loop over layers and partitions to access density
        in.getPartitionNames(partitions);
        const string dstr = "density";
        fieldName.push_back(dstr);
        
        BOOST_FOREACH (const string &part, partitions) {
            vertexFields.push_back(part);
        }
        BOOST_FOREACH (const string &partition, partitions) {

            if (!matchString(partition, fieldName)) {
                continue;
            }

            in.getScalarLayerNames(layers, partition);

            BOOST_FOREACH (const string &scalarLayer, layers) {

                if (!matchString(scalarLayer, fieldName)) {
                    continue;
                }  
                // if our f3d file contains float data
                Field<float>::Vec fScalarFields;
                fScalarFields = in.readScalarLayers<float>(partition, scalarLayer);

                if (fScalarFields.size() > 0) {

                    ReadLayersAndSetupBounds<float>(fScalarFields, minBound, maxBound, mapping); 
                }

                // if our f3d file contains double data
                Field<double>::Vec dScalarFields;
                dScalarFields = in.readScalarLayers<double>(partition, scalarLayer);

                if (dScalarFields.size() > 0) {

                    ReadLayersAndSetupBounds<double>(dScalarFields, minBound, maxBound, mapping); 

                }

                // if our f3d file contains double data
                Field<half>::Vec hScalarFields;
                hScalarFields = in.readScalarLayers<half>(partition, scalarLayer);

                if (hScalarFields.size() > 0) {

                    ReadLayersAndSetupBounds<half>(hScalarFields, minBound, maxBound, mapping); 

                }

            } 
        } 
        if (doblur != 0) {
            cout << "MotionBegin [" << shutter_open << " " << shutter_close << "]" << endl;
        }
        cout << "Volume \"blobbydso:F3DImplicitField.so\" [" 
             << minBound.x << " " 
             << maxBound.x << " " 
             << minBound.y << " " 
             << maxBound.y << " " 
             << minBound.z << " " 
             << maxBound.z << "]"
             << "[2 2 2] \"constant float[2] blobbydso:floatargs\" [" 
             << blur << " " << bbox_mod << " " << blur_cubic << " " << field_cubic << "]" 
             << "  \"constant string[1] blobbydso:stringargs\" [ \"" 
             << filename << "\" ]"
             << " \"constant float blobbydso:threshold\" ["
             << threshold << "] ";

        vector<string>::const_iterator vit = vertexFields.begin();
        for (; vit != vertexFields.end(); ++vit)
        {
            if (*vit != "vel" && *vit != "density") {
                cout << "\"varying float " << *vit << "\" [0 0 0 0 0 0 0 0]";
            }
        }
        cout << endl;

        if (doblur != 0) {
            cout << "Volume \"blobbydso:F3DImplicitField.so\" [" 
                 << minBound.x - blur_distance << " " 
                 << maxBound.x + blur_distance << " " 
                 << minBound.y - blur_distance << " " 
                 << maxBound.y + blur_distance << " " 
                 << minBound.z - blur_distance << " " 
                 << maxBound.z + blur_distance << "]"
                 << "[2 2 2] \"constant float[2] blobbydso:floatargs\" [" 
                 << blur << " " << bbox_mod << " " << blur_cubic << " " << field_cubic << "]" 
                 << "  \"constant string[1] blobbydso:stringargs\" [ \"" 
                 << filename << "\" ]"
                 << " \"constant float blobbydso:threshold\" ["
                 << threshold << "] ";

            vector<string>::const_iterator vit = vertexFields.begin();
            for (; vit != vertexFields.end(); ++vit)
            {
                if (*vit != "vel" && *vit != "density") {
                    cout << "\"varying float " << *vit << "\" [0 0 0 0 0 0 0 0]";
                }
            }
            cout << endl;
            cout << "MotionEnd" << endl;

        }

        printf("\377");
        fflush(stdout);
    }

    printf("\377");
    fflush(stdout);
    return 0;
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

template <typename T>
void ReadLayersAndSetupBounds(typename Field<T>::Vec sFields, 
                              V3d &minBound,
                              V3d &maxBound,
                              FieldMapping *mapping)
{
    V3d minvec;
    V3d maxvec;

    BOOST_FOREACH (typename Field<T>::Ptr field, sFields) {

        /* Get info from cache */
        Box3i field_extents = field->extents();
        const V3i minres = field_extents.min;
        const V3i maxres = field_extents.max;
        minres.getValue(minvec.x, minvec.y, minvec.z);
        maxres.getValue(maxvec.x, maxvec.y, maxvec.z);
    
        mapping = field->mapping().get();

        // derive object-space bounds from extents
        V3d big((V3d)maxres + V3d(1));
        maxBound = xformPoint(big, mapping);

        V3d small((V3d)minres + V3d(1));
        minBound = xformPoint(small, mapping);

    }
}

V3d xformPoint(V3d &Pt, FieldMapping *mapping)
{
    M44d vTw;
    MatrixFieldMapping::Ptr xmapping = 
        boost::dynamic_pointer_cast<MatrixFieldMapping>(mapping);
    vTw = xmapping->voxelToWorld();

    M44d big;
    big.setTranslation(Pt);
    M44d mres;
    mres.multiply(big, vTw, mres);
    V3d ret = mres.translation();

    return ret;
}


