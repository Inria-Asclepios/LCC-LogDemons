#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <string>

#include <tclap/CmdLine.h>
#include <tinyxml.h>

#include <itkTransform.h>
#include <itkDisplacementFieldTransform.h>
#include <itkStationaryVelocityFieldTransform.h>
#include <itkGeneralTransform.h>
#include <itkImageIOBase.h>
#include <itkImageIOFactory.h>
#include <itkImageFileReader.h>

#include <itkTransformToDeformationFieldSource.h>

#include <mipsInrimageImageIOFactory.h>

#include "rpiCommonTools.hxx"

/**
 * Fuse a list of transformations into a single transformation. The list of transformations must
 * be stored into a valid XML file similar to this example:
 *
 *   <?xml version="1.0" encoding="UTF-8"?>
 *   <listoftransformations>
 *       <transformation>
 *           <type>linear</type>
 *           <path>t0.txt</path>
 *       </transformation>
 *       <transformation>
 *           <type>displacementfield</type>
 *           <path>t1.nii.gz</path>
 *       </transformation>
 *   </listoftransformations>
 *
 * In this example, we have 2 transformations. The first one (denoted T0) is a linear transformation
 * stored into the text file "t0.txt". The second transformation (denoted T1) is a displacement
 * field stored into the -- image -- file "t1.nii.gz". The computed transformation is given by
 * T0 o T1. Thus, the order the transformations written into the XML file is important since the
 * composition of two transformations is usually not commutative. The XML file can contain any
 * strictly positive number of transformations.
 *
 * The meaning of the different tags is:
 *
 *  -Tag <listoftransformations> contains a set of transformations. This tag is mandatory.
 *
 *  -Tag <transformation> contains a transformation. This tag must be an element of the tag
 *   <listoftransformations>. This tag is mandatory.
 *
 *  -Tag <type> contains the type of the considered transformation. The available types are "linear",
 *   "displacemendfield", and "stiaionaryvelocityfield". This tag must be an element of the tag
 *   <transformation>. This tag is mandatory.
 *
 *  -Tag <path> contains the path to the file containing the transformation. A displacement field
 *   or a stationary velocity field must be stored into an image file (e.g. Nifty, Analyse, etc.).
 *   This tag must be an element of the tag <transformation>. This tag is mandatory.
 *
 *  -Tag <invert> indicates if the considered transformation has to be inverted before the fusion.
 *   The value "1" induce the transformation inversion ; the value "0" means no inversion. The only
 *   transformations that can be inverted yet are the linear transformations. This tag is optional.
 *   By default, no inversion is performed.
 *
 * The type of the output transformation depends on the transformations of the list. A list of
 * linear transformation will be fused into a linear transformation. If the list contains at least
 * one displacement field, the output transformation will be a displacement field. If the option
 * "--force-displacement-field" is used, the output transformation will be a displacement field.
 * One should choose the output file extension according to the transformation computed. Indeed, if
 * the output transformation is a displacement field, only image format are supported (e.g. .nii).
 * No checking on the file extension will be performed since the output transformation type is
 * easily predictable.
 *
 * @author Vincent Garcia
 * @date   2011/05/31
 */


/**
 * Transformation type
 */
enum Transformation {
    LINEAR,                    /** Linear                    */
    DISPLACEMENT_FIELD,        /** Displacement field        */
    STATIONARY_VELOCITY_FIELD  /** Stationary velocity field */
};


/**
 * Structure containing the IO parameters, the verbosity,
 * and eventually the output transformation type.
 */
struct Param
{
    std::string inputFile;
    std::string outputFile;
    std::string geometryFile;
    bool        forceDisplacementField;
    bool        verbose;
};


/**
 * Structure containing the information of the transformations (type, path, inversion).
 */
struct Data
{
    std::vector<Transformation> type;
    std::vector<std::string>    path;
    std::vector<bool>           invert;
};


/**
 * Returns the transformation type as string.
 * @param  transformation transformation
 * @return string describing the transformation
 */
std::string getTransformationAsString(Transformation transformation)
{
    if (transformation==LINEAR)
        return "linear";
    else if (transformation==DISPLACEMENT_FIELD)
        return "displacement field";
    else if (transformation==STATIONARY_VELOCITY_FIELD)
        return "stationary velocity field";
    else
        return "transformation not recognized";
}


/**
 * Parses the command line arguments and deduces the corresponding Param structure.
 * @param  argc   number of arguments
 * @param  argv   array containing the arguments
 * @param  param  structure of parameters
 */
void parseParameters(int argc, char** argv, struct Param & param)
{

    // Program description
    std::string description = "\b\b\bDESCRIPTION\n";

    description += "Fuse a list of transformations into a single transformation. The list of transformations must ";
    description += "be stored into a valid XML file similar to this example:\n";

    description += "  <?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    description += "  <listoftransformations>\n";
    description += "      <transformation>\n";
    description += "          <type>linear</type>\n";
    description += "          <path>t0.txt</path>\n";
    description += "          <invert>0</invert>\n";
    description += "      </transformation>\n";
    description += "      <transformation>\n";
    description += "          <type>displacementfield</type>\n";
    description += "          <path>t1.nii.gz</path>\n";
    description += "          <invert>0</invert>\n";
    description += "      </transformation>\n";
    description += "  </listoftransformations>\n";

    description += "In this example, we have 2 transformations. The first one (denoted T0) is a linear transformation ";
    description += "stored into the text file \"t0.txt\". The second transformation (denoted T1) is a displacement ";
    description += "field stored into the -- image -- file \"t1.nii.gz\". The computed transformation is given by ";
    description += "T0 o T1. Thus, the order the transformations written into the XML file is important since the ";
    description += "composition of two transformations is usually not commutative. The XML file can contain any ";
    description += "strictely positive number of transformations.\n";

    description += "The meaning of the different tags is:\n";

    description += " -Tag <listoftransformations> contains a set of transformations. This tag is mandatory.\n";

    description += " -Tag <transformation> contains a transformation. This tag must be an element of the tag ";
    description += "  <listoftransformations>. This tag is mandatory.\n";

    description += " -Tag <type> contains the type of the considered transformation. The available types are \"linear\", ";
    description += "  \"displacemendfield\", and \"stiaionaryvelocityfield\". This tag must be an element of the tag ";
    description += "  <transformation>. This tag is mandatory.\n";

    description += " -Tag <path> contains the path to the file containing the transformation. A displacement field ";
    description += "  or a stationary velocity field must be stored into an image file (e.g. Nifty, Analyse, etc.). ";
    description += "  This tag must be an element of the tag <transformation>. This tag is mandatory.\n";

    description += " -Tag <invert> indicates if the considered transformation has to be inverted before the fusion. ";
    description += "  The value \"1\" induces the transformation inversion ; the value \"0\" means no inversion. ";
    description += "  This tag is optional. By default, no inversion is performed.\n";

    description += "The type of the output transformation depends on the transformations of the list. A list of ";
    description += "linear transformations will be fused into a linear transformation. If the list contains at least ";
    description += "one displacement field, the output transformation will be a displacement field. If the option ";
    description += "\"--force-displacement-field\" is used, the ouput transformation will be a displacement field. ";
    description += "One should choose the output file extension according to the transformation computed. Indeed, if ";
    description += "the output transformation is a displacement field, only image format are supported (e.g. .nii). ";
    description += "No checking on the file extension will be performed since the output transformation type is ";
    description += "easily predictable.";

    description += "";
    description += "\nAuthors : Vincent Garcia";

    // Option description
    std::string dInput     = "Path to the list of transformations (XML file).";

    std::string dOutput    = "Path to the output transformation file.";

    std::string dGeometry  = "Path to the image containing the geometry of the output transformation. ";
    dGeometry             += "The geometry is used only if the transformation computed is either a ";
    dGeometry             += "displacement field or a stationary velocity field.";

    std::string dForce     = "Force the output transformation to be a displacement field.";

    std::string dVerbose   = "Verbose mode.";

    try {

        // Define the command line parser
        TCLAP::CmdLine cmd( description, ' ', "1.0", true);

        // Set options
        TCLAP::SwitchArg              aVerbose(  "",  "verbose",                  dVerbose, cmd, false);
        TCLAP::SwitchArg              aForce(    "f", "force-displacement-field", dForce,   cmd, false);
        TCLAP::ValueArg<std::string>  aGeometry( "g", "geometry",         dGeometry, false, "", "string", cmd );
        TCLAP::ValueArg<std::string>  aOutput(   "o", "output-transform", dOutput,   true,  "", "string", cmd );
        TCLAP::ValueArg<std::string>  aInput(    "i", "input-list",       dInput,    true,  "", "string", cmd );

        // Parse the command line
        cmd.parse( argc, argv );

        // Set the parameters
        param.inputFile              = aInput.getValue();
        param.outputFile             = aOutput.getValue();
        param.geometryFile           = aGeometry.getValue();
        param.forceDisplacementField = aForce.getValue();
        param.verbose                = aVerbose.getValue();
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
        throw std::runtime_error("Unable to parse the command line arguments.");
    }
}


/**
 * Parses the XML file anf fill the data structure.
 * @param fileName  path to the XML file
 * @param data      Data structure
 */
void parseXML(const char* fileName, struct Data & data)
{

    // Load XML document
    TiXmlDocument doc(fileName);
    bool loadOkay = doc.LoadFile();
    if (loadOkay)
    {

        // Variables
        TiXmlElement *pRoot, *pParm, *pType, *pPath, *pInvert;

        // Get XML root
        pRoot = doc.FirstChildElement( "listoftransformations" );
        if ( pRoot )
        {
            // Get first transformation
            pParm = pRoot->FirstChildElement("transformation");
            while ( pParm )
            {
                // Get transformation type, transformation path, and the inversion
                pType   = pParm->FirstChildElement("type");
                pPath   = pParm->FirstChildElement("path");
                pInvert = pParm->FirstChildElement("invert");

                // Process transformation
                if ( pType && pPath )
                {

                    // Type and path
                    std::string type(pType->GetText());
                    std::string path(pPath->GetText());

                    // Process type
                    if (type.compare("linear")==0)
                        data.type.push_back(LINEAR);
                    else if (type.compare("displacementfield")==0)
                        data.type.push_back(DISPLACEMENT_FIELD);
                    else if (type.compare("stationaryvelocityfield")==0)
                        data.type.push_back(STATIONARY_VELOCITY_FIELD);
                    else
                    {
                        std::string warning;
                        warning += "Warning : the transformation type \"" + type + "\" is not ";
                        warning += "supported. The transformation " + path + " will be ignored.";
                        std::cout << warning << std::endl;
                        pParm = pParm->NextSiblingElement("transformation");
                        break;
                    }

                    // Process path
                    data.path.push_back(path);

                    // Process inversion
                    if ( pInvert )
                    {
                        std::string invert(pInvert->GetText());
                        if (invert.compare("1")==0)
                        {
                            if ( type.compare("stationaryvelocityfield")==0 )
                            {
                                std::string warning;
                                warning += "Warning : a transformation of type \"" + type + "\" cannot ";
                                warning += "be inverted. The inversion will not be performed.";
                                std::cout << warning << std::endl;
                                data.invert.push_back(false);
                            }
                            else
                                data.invert.push_back(true);
                        }
                        else
                            data.invert.push_back(false);
                    }
                    else
                        data.invert.push_back(false);
                }
                else
                    throw std::runtime_error("Transformation must contain tags \"type\" and \"path\".");

                // Get next transformation
                pParm = pParm->NextSiblingElement("transformation");
            }
        }
    }
    else
        throw std::runtime_error( std::string("Failed to load file ") + fileName + "." );
}


/**
 * Prints the Param structure.
 * @param param    Param structure
 * @param verbose  prints if and only if verbose is true
 */
void printParam(struct Param & param, bool verbose)
{
    if (!verbose)
        return;

    std::cout << std::endl;
    std::cout << "PARAMETERS" << std::endl << std::endl;
    std::cout << "  Input list               : " << param.inputFile << std::endl;
    std::cout << "  Output transform         : " << param.outputFile << std::endl;
    if (param.geometryFile.compare("")!=0)
        std::cout << "  Geometry                 : " << param.geometryFile << std::endl;
    if (param.forceDisplacementField)
        std::cout << "  Force displacement field : true" << std::endl;
    else
        std::cout << "  Force displacement field : false" << std::endl << std::endl;
}


/**
 * Prints the Data structure.
 * @param data     Data structure
 * @param verbose  prints if and only if verbose is true
 */
void printData(struct Data & data, bool verbose)
{
    if (!verbose)
        return;

    std::cout << std::endl;
    std::cout << "TRANSFORMATION LIST" << std::endl << std::endl;
    for (unsigned int i=0; i<data.type.size(); i++)
    {
        std::cout << "  Transformation " << i << std::endl;

        if (data.type[i]==LINEAR)
            std::cout << "  Type   : linear" << std::endl;
        else if (data.type[i]==DISPLACEMENT_FIELD)
            std::cout << "  Type   : displacement field" << std::endl;
        else if (data.type[i]==STATIONARY_VELOCITY_FIELD)
            std::cout << "  Type   : stationary velocity field" << std::endl;

        std::cout << "  Path   : " << data.path[i]   << std::endl;

        std::cout << "  Invert : " << data.invert[i] << std::endl << std::endl;
    }
}


/**
 * Build the list of transformations.
 * @param  data data structure
 * @return list of transformations
 */
template<class TScalarType>
typename itk::GeneralTransform<TScalarType,3>::Pointer
buildListOfTransformations(struct Data & data)
{
    // Type definition
    typedef itk::GeneralTransform<TScalarType,3>               ListType;
    typedef itk::Transform<TScalarType,3,3>                    LinearType;
    typedef itk::MatrixOffsetTransformBase<TScalarType,3,3>    MatrixType;
    typedef itk::DisplacementFieldTransform<TScalarType>       DFType;
    typedef itk::StationaryVelocityFieldTransform<TScalarType> SVFType;

    // Create and fill the list of transformations
    typename ListType::Pointer list = ListType::New();
    for (unsigned int i=0; i<data.type.size(); i++)
    {

        if (data.type[i]==LINEAR)
        {

            typename LinearType::Pointer linear = rpi::readLinearTransformation<TScalarType>(data.path[i]);

            if (data.invert[i]==false)
                list->InsertTransform( linear.GetPointer() );
            else
            {
                MatrixType * matrix = dynamic_cast<MatrixType *>(linear.GetPointer());
                if (matrix==0)
                    throw std::runtime_error("Cannot cast the transformation into an itk::MatrixOffsetTransformBase.");
                typename MatrixType::Pointer inverse = MatrixType::New();
                inverse->SetCenter(matrix->GetCenter());
                matrix->GetInverse(inverse);
                list->InsertTransform( inverse.GetPointer() );
            }

        }
        else if (data.type[i]==DISPLACEMENT_FIELD)
        {
            typename DFType::Pointer field = rpi::readDisplacementField<TScalarType>(data.path[i]);

            if (data.invert[i]==true)
                field->GetInverse(field);
            list->InsertTransform( field.GetPointer() );


        }
        else if (data.type[i]==STATIONARY_VELOCITY_FIELD)
        {
            typename SVFType::Pointer field = rpi::readStationaryVelocityField<TScalarType>(data.path[i]);

            if (data.invert[i]==true)
                field->GetInverse(field);
            list->InsertTransform( field.GetPointer() );
        }
        else
            throw std::runtime_error("Transformation not supported.");
    }
    return list;
}


/**
 * Decides the type of transformation computed.
 * @param  list list of transformation
 * @return transformation type as enum
 */
template<class TScalarType>
Transformation decideOuputTransformationType(typename itk::GeneralTransform<TScalarType,3>::Pointer list)
{
    if (list->IsLinear())
        return LINEAR;
    else
        return DISPLACEMENT_FIELD;
}


/**
 * Generates linear transformation from the transformation list and exports it.
 * @param list      list of transformations
 * @param fileName  output file name
 */
template<class TScalarType>
void exportLinearTransformation(itk::GeneralTransform<TScalarType,3> * list, std::string fileName)
{
    // Type definition
    typedef itk::MatrixOffsetTransformBase<TScalarType, 3, 3>  MatrixOffsetTransformType;
    typedef itk::AffineTransform<TScalarType,3>                AffineTransformType;

    // Get the MatrixOffsetTransformBase from the list
    typename MatrixOffsetTransformType::Pointer matrix = list->GetGlobalLinearTransform();

    // Create the AffineTransform from the MatrixOffsetTransformBase
    typename AffineTransformType::Pointer affine = AffineTransformType::New();
    affine->SetMatrix(matrix->GetMatrix());
    affine->SetOffset(matrix->GetOffset());
    affine->SetCenter(matrix->GetCenter());

    // Write the AffineTransform
    rpi::writeLinearTransformation<TScalarType,3>(affine, fileName);
}


/**
 * Generates a displacement field from the transformation list and exports it as image.
 * @param list              list of transformations
 * @param data              data structure
 * @param geometryFileName  path to the image containing the geometry
 * @param fileName          output file name
 */
template<class TScalarType>
void exportDFTransformation(itk::GeneralTransform<TScalarType,3> * list, struct Data data, std::string geometryFileName, std::string fileName)
{

    // Type definition
    typedef  itk::GeneralTransform<TScalarType,3>                                    TransformListType;
    typedef  itk::DisplacementFieldTransform< TScalarType, 3 >                       DFType;
    typedef  itk::StationaryVelocityFieldTransform< TScalarType, 3 >                 SVFType;
    typedef  typename DFType::VectorFieldType                                        VectorFieldType;
    typedef  itk::TransformToDeformationFieldSource< VectorFieldType, TScalarType >  GeneratorType;

    // Create a field generator
    typename GeneratorType::Pointer fieldGenerator = GeneratorType::New();
    fieldGenerator->SetTransform( list );

    // Sets the geometry of the displacement field
    if (geometryFileName.compare("")!=0)
    {
        // Type definition
        typedef  itk::ImageIOBase                     ImageIOBase;
        typedef typename GeneratorType::PointType     OriginType;
        typedef typename GeneratorType::SpacingType   SpacingType;
        typedef typename GeneratorType::DirectionType DirectionType;
        typedef typename GeneratorType::RegionType    RegionType;

        // Read image header
        typename ImageIOBase::Pointer imageBase = rpi::readImageInformation( geometryFileName );

        // Check the image dimension
        if (imageBase->GetNumberOfDimensions()!=3)
            throw std::runtime_error("Only 3D images are supported yet.");

        // Create and fill origin, spacing, direction, and region from image header
        OriginType    origin;
        SpacingType   spacing;
        DirectionType direction;
        RegionType    region;
        for (int i=0; i<3; i++)
        {
            origin[i]  = imageBase->GetOrigin(i);
            spacing[i] = imageBase->GetSpacing(i);
            region.SetSize(i,imageBase->GetDimensions(i));
            for (int j=0; j<3; j++)
                direction(i,j) = imageBase->GetDirection(j)[i];
        }

        // Set the geometry of the field generator
        fieldGenerator->SetOutputOrigin(    origin );
        fieldGenerator->SetOutputSpacing(   spacing );
        fieldGenerator->SetOutputDirection( direction );
        fieldGenerator->SetOutputRegion(    region );
    }
    else
    {
        // Locate the first non linear transformation in the list to get its geometry
        unsigned int index;
        for (index=0; index<data.type.size(); index++)
            if ( data.type[index]==DISPLACEMENT_FIELD || data.type[index]==STATIONARY_VELOCITY_FIELD )
                break;

        // Stop if there is no non linear transformation in the list
        if (index==data.type.size())
            throw std::runtime_error("No geometry was found in the transformation list.");

        // Get the geometry for the field located
        const DFType  * df  = dynamic_cast<const DFType*>( list->GetTransform(index).GetPointer());
        const SVFType * svf = dynamic_cast<const SVFType*>(list->GetTransform(index).GetPointer());
        if (df!=0)
        {
            typename VectorFieldType::ConstPointer container = df->GetParametersAsVectorField();
            fieldGenerator->SetOutputOrigin(    container->GetOrigin() );
            fieldGenerator->SetOutputSpacing(   container->GetSpacing() );
            fieldGenerator->SetOutputDirection( container->GetDirection() );
            fieldGenerator->SetOutputRegion(    container->GetLargestPossibleRegion() );
        }
        else if (svf!=0)
        {
            typename VectorFieldType::ConstPointer container = svf->GetParametersAsVectorField();
            fieldGenerator->SetOutputOrigin(    container->GetOrigin() );
            fieldGenerator->SetOutputSpacing(   container->GetSpacing() );
            fieldGenerator->SetOutputDirection( container->GetDirection() );
            fieldGenerator->SetOutputRegion(    container->GetLargestPossibleRegion() );
        }
        else
            throw std::runtime_error("No geometry can be extracted from the transformation.");
    }

    // Update the field generator
    try
    {
        fieldGenerator->Update();
    }
    catch( itk::ExceptionObject& err )
    {
        throw std::runtime_error( "Could not generate a displacement field from a linear transformation." );
    }

    // Gets the displacement field (image)
    typename VectorFieldType::Pointer container = fieldGenerator->GetOutput();
    container->DisconnectPipeline();

    // Create the displacement field (transformation)
    typename DFType::Pointer field = DFType::New();
    field->SetParametersAsVectorField( static_cast<typename VectorFieldType::ConstPointer>( container.GetPointer() ) );

    // Write displacement field
    rpi::writeDisplacementFieldTransformation<TScalarType,3>(field, fileName );
}


/**
 * Main function.
 */
int main(int argc, char** argv)
{

    // Type definition
    typedef float                               ScalarType;
    typedef itk::GeneralTransform<ScalarType,3> TransformListType;

    // Structures
    struct Param param;
    struct Data  data;

    // Allows the executable to read and write Inrimage
    itk::InrimageImageIOFactory::RegisterOneFactory();

    try{

        // Parse command line options
        parseParameters(argc, argv, param);
        printParam(param, param.verbose);

        // Parse XML file
        parseXML(param.inputFile.c_str(), data);
        printData(data, param.verbose);

        // Build the list of transformations
        TransformListType::Pointer list = buildListOfTransformations<ScalarType>(data);

        // Print information
        if (param.verbose)
            std::cout << "INFORMATIONS" << std::endl;

        // Decide the type of the output transformation
        Transformation output;
        if (param.forceDisplacementField)
            output = DISPLACEMENT_FIELD;
        else
            output = decideOuputTransformationType<ScalarType>(list);
        if (param.verbose)
            std::cout << "  Output transformation type : " << getTransformationAsString(output) << std::endl;

        // Switch on the output transformation type
        if (param.verbose)
            std::cout << "  Fusion progress            : " << std::flush;
        if (output==LINEAR)
            exportLinearTransformation<ScalarType>(list, param.outputFile);
        else if (output==DISPLACEMENT_FIELD)
            exportDFTransformation<ScalarType>(list, data, param.geometryFile, param.outputFile);
        else
            throw std::runtime_error( "Transformation type not supported yet." );
        if (param.verbose)
            std::cout << "done" << std::endl;

    }
    catch( std::exception& e )
    {
        std::cerr << std::endl << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    };

    return EXIT_SUCCESS;
}
