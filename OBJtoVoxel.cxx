// #include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkAppendFilter.h>
#include <vtkPoints.h>
#include <vtkVoxel.h>
// #include <vtkNamedColors.h>
#include <vtkDataSetMapper.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkCellData.h>
#include <vtkOBJReader.h>
#include <vtkImageData.h>
#include <vtkPolyDataMapper.h>
#include <vtkUnstructuredGrid.h>
// #include <vtkUnstructuredGridVolumeMapper.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkProperty.h>
#include <vtkVoxelModeller.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkMultiThreshold.h>
#include <vtkImageDataToPointSet.h>
#include <vtkImageShrink3D.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
// #include <vtkRenderWindow.h>
// #include <vtkRenderWindowInteractor.h>
// #include <vtkRenderer.h>

#include <string>
#include <iostream>

int main(int argc, char* argv[])
{
  // Parse command line arguments
  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << "Filename(.obj) e.g trumpet.obj"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];
  std::cout << "Read File : " << filename << std::endl;
  vtkNew<vtkOBJReader> reader;
  reader->SetFileName(filename.c_str());
  reader->Update();
  std::cout << "Complete to read file : " << std::endl;
   // Write the file
  vtkNew<vtkXMLPolyDataWriter> PolyWriter;
  PolyWriter->SetFileName("obj.vtp");
  PolyWriter->SetInputData(reader->GetOutput());
  // Optional - set the mode. The default is binary.
  // writer->SetDataModeToBinary();
  // writer->SetDataModeToAscii();
  PolyWriter->Write();

  vtkNew<vtkSphereSource> sphereSource;
  sphereSource->SetRadius(20);
  sphereSource->SetPhiResolution(30);
  sphereSource->SetThetaResolution(30);
  auto pd = sphereSource->GetOutput();
  sphereSource->Update();

  vtkNew<vtkPolyData> polydata;

  // vtkNew<vtkImageData> whiteImage;

  // std::cout << "GetBounds" << std::endl;
  double bounds[6];
  reader->GetOutput()->GetBounds(bounds);
  std::cout << "Bounds[0,1,2,3,4,5] : (" << bounds[0] << "," << bounds[1] << "," << bounds[2] << "," << bounds[3] << "," << bounds[4] << "," << bounds[5] << ")" << std::endl;
  bounds[0] = bounds[0] * 0.2;
  bounds[1] = bounds[1] * 0.7;
  bounds[2] = bounds[2] * 0.7;
  bounds[3] = bounds[3] * 0.6;
  bounds[4] = -5.0;
  bounds[5] = 40;
  // // pd->GetBounds(bounds);
  double spacing[3]; // desired volume spacing
  spacing[0] = 2.0;
  spacing[1] = 2.0;
  spacing[2] = 2.0;
  // whiteImage->SetSpacing(spacing);

  // compute dimensions  
  std::cout << "Compute Dimensions" << std::endl;
  int dim[3];
  for (int i = 0; i < 3; i++)
  {
    dim[i] = static_cast<int>(
        ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
  }
  std::cout << "Dimensions[x,y,z] : (" << dim[0] << "," << dim[1] << "," << dim[2] << ")" << std::endl;
  // whiteImage->SetDimensions(dim);
  // whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);
  // double origin[3];
  // origin[0] = bounds[0] + spacing[0] / 2;
  // origin[1] = bounds[2] + spacing[1] / 2;
  // origin[2] = bounds[4] + spacing[2] / 2;
  // whiteImage->SetOrigin(origin);
  // whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
  // // Fill the image with foreground voxels:
  // std::cout << "Fill the image with foreground voxels" << std::endl;
  // unsigned char inval = 255;
  // unsigned char outval = 0;
  // vtkIdType count = whiteImage->GetNumberOfPoints();
  // for (vtkIdType i = 0; i < count; ++i)
  // {
  //   whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
  // }

  // // polygonal data --> image stencil:
  // std::cout << "Polygonal Data --> Image Stencil" << std::endl;
  // vtkNew<vtkPolyDataToImageStencil> pol2stenc;
  // pol2stenc->SetInputData(reader->GetOutput());
  // // pol2stenc->SetInputData(pd);
  // pol2stenc->SetOutputOrigin(origin);
  // pol2stenc->SetOutputSpacing(spacing);
  // pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
  // std::cout << "Update : Polygonal Data --> Image Stencil" << std::endl;
  // pol2stenc->Update();
  
  // // Cut the corresponding white image and set the background:
  // std::cout << "Cut the corresponding white image and set the background" << std::endl;
  // vtkNew<vtkImageStencil> imgstenc;
  // imgstenc->SetInputData(whiteImage);
  // imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
  // imgstenc->ReverseStencilOff();
  // imgstenc->SetBackgroundValue(outval);
  // imgstenc->Update();

  // // Output
  // std::cout << "Output ImageData" << std::endl;
  // vtkNew<vtkMetaImageWriter> writer;
  // writer->SetFileName("Volume.mhd");
  // writer->SetInputData(imgstenc->GetOutput());
  // writer->Write();

  std::cout << "Set VoxelModeller" << std::endl;
  vtkNew<vtkVoxelModeller> voxelModeller;
  voxelModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
  voxelModeller->SetModelBounds(bounds);
  voxelModeller->SetScalarTypeToFloat();
  voxelModeller->SetMaximumDistance(0.1);
  // voxelModeller->SetForegroundValue(10.0);
  // voxelModeller->SetBackgroundValue(0.0);
  voxelModeller->SetInputConnection(reader->GetOutputPort());
  voxelModeller->Update();
  
  vtkNew<vtkImageData> volume;
  volume->DeepCopy(voxelModeller->GetOutput());
  // Output
  std::cout << "Output VoxelModeller ImageData" << std::endl;
  vtkNew<vtkMetaImageWriter> writerImage;
  writerImage->SetFileName("VoxelModeller.mhd");
  writerImage->SetRAWFileName("VoxelModeller.raw");
  writerImage->SetInputData(voxelModeller->GetOutput());
  writerImage->Write();

  // vtkNew<vtkImageShrink3D> shrink;
  // shrink->SetShrinkFactors(4, 4, 4);
  // shrink->SetInputConnection(imgstenc->GetOutputPort());
  // shrink->Update();

  std::cout << "Convert ImageData to PointSet" << std::endl;
  vtkNew<vtkImageDataToPointSet> imageDataToPointSet;
  // imageDataToPointSet->SetInputConnection(imgstenc->GetOutputPort());
  // imageDataToPointSet->SetInputConnection(shrink->GetOutputPort());
  imageDataToPointSet->SetInputData(voxelModeller->GetOutput());
  imageDataToPointSet->Update();

  // Extract voxels on the border between the inside and outside.
  vtkNew<vtkMultiThreshold> threshold;
  // Inside points have one or more points above the isosurface
  int insideId = threshold->AddIntervalSet(
      -0.01, 0.01, vtkMultiThreshold::CLOSED, vtkMultiThreshold::CLOSED,
      vtkDataObject::FIELD_ASSOCIATION_POINTS, "ImageScalars", 0, 0);
  // Border points have points that straddle the boundary
  int borderId = threshold->AddIntervalSet(
      0.99, 1.01, vtkMultiThreshold::OPEN, vtkMultiThreshold::OPEN,
      vtkDataObject::FIELD_ASSOCIATION_POINTS, "ImageScalars", 0, 0);

  threshold->SetInputData(imageDataToPointSet->GetOutput());

  // Select the intervals to be output
  threshold->OutputSet(insideId);
  threshold->OutputSet(borderId);
  threshold->Update();
  
  // vtkNew<vtkSelectEnclosedPoints> select;
  // select->SetInputData(reader->GetOutput());
  // select->SetSurfaceData(reader->GetOutput());

  // vtkNew<vtkMultiThreshold> threshold;
  // // Outside points have a 0 value in ALL points of a cell
  // int outsideId = threshold->AddBandpassIntervalSet(
  //     0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "SelectedPoints", 0, 1);
  // // Inside points have a 1 value in ALL points of a cell
  // int insideId = threshold->AddBandpassIntervalSet(
  //     1, 1, vtkDataObject::FIELD_ASSOCIATION_POINTS, "SelectedPoints", 0, 1);
  // // Border points have a 0 or a 1 in at least one point of a cell
  // int borderId = threshold->AddIntervalSet(
  //     0, 1, vtkMultiThreshold::OPEN, vtkMultiThreshold::OPEN,
  //     vtkDataObject::FIELD_ASSOCIATION_POINTS, "SelectedPoints", 0, 0);

  // threshold->SetInputConnection(select->GetOutputPort());

  // // Select the intervals to be output
  // threshold->OutputSet(outsideId);
  // threshold->OutputSet(insideId);
  // threshold->OutputSet(borderId);
  // threshold->Update();


  vtkNew<vtkAppendFilter> appendFilter;
  appendFilter->AddInputData(reader->GetOutput());
  appendFilter->Update();

  
  vtkNew<vtkUnstructuredGrid> unstructuredGrid_inside;
  // unstructuredGrid->ShallowCopy(appendFilter->GetOutput());
  unstructuredGrid_inside->ShallowCopy(dynamic_cast<vtkUnstructuredGrid*>(
      vtkMultiBlockDataSet::SafeDownCast(
          threshold->GetOutput()->GetBlock(insideId))
          ->GetBlock(0)));

  vtkNew<vtkUnstructuredGrid> unstructuredGrid_voxel_inside;
  vtkNew<vtkVoxel> voxel_inside;
  vtkNew<vtkPoints> points_inside;
  vtkNew<vtkFloatArray> FloatArray_inside;
  std::cout << "Convert Hexahedron to Voxel (Inside Voxel)" << std::endl;
  int idCount = 0;
  for (int cellIdx = 0; cellIdx < unstructuredGrid_inside->GetNumberOfCells(); ++cellIdx){
    float objValue = 0.0;
    FloatArray_inside->InsertNextValue(objValue);
    vtkNew<vtkIdList> idList;
    vtkNew<vtkIdList> idListTemp;
    unstructuredGrid_inside->GetCellPoints(cellIdx,idList);
    // std::cout << "idList numbers : " << idList->GetNumberOfIds() << std::endl;
    double pointsTemp[8][3];
    for(int pointIdx = 0; pointIdx < idList->GetNumberOfIds();pointIdx++){
      // std::cout << "idList [" << pointIdx << "] : " << idList->GetId(pointIdx) << std::endl;
      unstructuredGrid_inside->GetPoint(idList->GetId(pointIdx),pointsTemp[pointIdx]);
      voxel_inside->GetPointIds()->SetId(pointIdx, idCount);
      idCount++;
      // std::cout << "Point of " << idList->GetId(pointIdx) << " : (" << pointsTemp[0] << "," << pointsTemp[1] << "," << pointsTemp[2] << ")" << std::endl;
    }
    points_inside->InsertNextPoint(pointsTemp[0]);
    points_inside->InsertNextPoint(pointsTemp[1]);
    points_inside->InsertNextPoint(pointsTemp[3]);
    points_inside->InsertNextPoint(pointsTemp[2]);
    points_inside->InsertNextPoint(pointsTemp[4]);
    points_inside->InsertNextPoint(pointsTemp[5]);
    points_inside->InsertNextPoint(pointsTemp[7]);
    points_inside->InsertNextPoint(pointsTemp[6]);
    unstructuredGrid_voxel_inside->InsertNextCell(voxel_inside->GetCellType(),voxel_inside->GetPointIds());
  }
  unstructuredGrid_voxel_inside->SetPoints(points_inside);
  unstructuredGrid_voxel_inside->GetCellData()->SetScalars(FloatArray_inside);

  std::cout << "Write Inside Voxel VTU file..." << std::endl;
  vtkNew<vtkXMLUnstructuredGridWriter> writerXML_inside;
  writerXML_inside->SetFileName("ObjToVoxel_inside.vtu");
  writerXML_inside->SetInputData(unstructuredGrid_voxel_inside);
  writerXML_inside->Write();


  vtkNew<vtkUnstructuredGrid> unstructuredGrid_border;
  // unstructuredGrid->ShallowCopy(appendFilter->GetOutput());
  unstructuredGrid_border->ShallowCopy(dynamic_cast<vtkUnstructuredGrid*>(
      vtkMultiBlockDataSet::SafeDownCast(
          threshold->GetOutput()->GetBlock(borderId))
          ->GetBlock(0)));

  vtkNew<vtkUnstructuredGrid> unstructuredGrid_voxel_border;
  vtkNew<vtkVoxel> voxel_border;
  vtkNew<vtkPoints> points_border;
  vtkNew<vtkFloatArray> FloatArray_border;
  std::cout << "Convert Hexahedron to Voxel (Border Voxel)" << std::endl;
  idCount = 0;
  for (int cellIdx = 0; cellIdx < unstructuredGrid_border->GetNumberOfCells(); ++cellIdx){
    float objValue = 0.0;
    FloatArray_border->InsertNextValue(objValue);
    vtkNew<vtkIdList> idList;
    vtkNew<vtkIdList> idListTemp;
    unstructuredGrid_border->GetCellPoints(cellIdx,idList);
    // std::cout << "idList numbers : " << idList->GetNumberOfIds() << std::endl;
    double pointsTemp[8][3];
    for(int pointIdx = 0; pointIdx < idList->GetNumberOfIds();pointIdx++){
      // std::cout << "idList [" << pointIdx << "] : " << idList->GetId(pointIdx) << std::endl;
      unstructuredGrid_border->GetPoint(idList->GetId(pointIdx),pointsTemp[pointIdx]);
      voxel_border->GetPointIds()->SetId(pointIdx, idCount);
      idCount++;
      // std::cout << "Point of " << idList->GetId(pointIdx) << " : (" << pointsTemp[0] << "," << pointsTemp[1] << "," << pointsTemp[2] << ")" << std::endl;
    }
    points_border->InsertNextPoint(pointsTemp[0]);
    points_border->InsertNextPoint(pointsTemp[1]);
    points_border->InsertNextPoint(pointsTemp[3]);
    points_border->InsertNextPoint(pointsTemp[2]);
    points_border->InsertNextPoint(pointsTemp[4]);
    points_border->InsertNextPoint(pointsTemp[5]);
    points_border->InsertNextPoint(pointsTemp[7]);
    points_border->InsertNextPoint(pointsTemp[6]);
    unstructuredGrid_voxel_border->InsertNextCell(voxel_border->GetCellType(),voxel_border->GetPointIds());
  }
  unstructuredGrid_voxel_border->SetPoints(points_border);
  unstructuredGrid_voxel_border->GetCellData()->SetScalars(FloatArray_border);

  std::cout << "Write Border Voxel VTU file..." << std::endl;
  vtkNew<vtkXMLUnstructuredGridWriter> writerXML_border;
  writerXML_border->SetFileName("ObjToVoxel_border.vtu");
  writerXML_border->SetInputData(unstructuredGrid_voxel_border);
  writerXML_border->Write();

  // vtkNew<vtkActor> actor;
  // actor->SetMapper(mapper);
  // actor->GetProperty()->SetDiffuseColor(actorColor.GetData());

  // vtkNew<vtkRenderer> renderer;
  // renderer->AddActor(actor);
  // renderer->SetBackground(backgroundColor.GetData());
  // renderer->ResetCamera();
  // renderer->GetActiveCamera()->Azimuth(30);
  // renderer->GetActiveCamera()->Elevation(30);
  // renderer->GetActiveCamera()->Dolly(1.5);
  // renderer->ResetCameraClippingRange();

  // vtkNew<vtkRenderWindow> renderWindow;
  // renderWindow->AddRenderer(renderer);
  // renderWindow->SetWindowName("ReadOBJ");

  // vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  // renderWindowInteractor->SetRenderWindow(renderWindow);

  // renderWindow->SetSize(640, 480);
  // renderWindow->Render();

  // renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
