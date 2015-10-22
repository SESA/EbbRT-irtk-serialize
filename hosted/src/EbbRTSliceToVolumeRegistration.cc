#include "EbbRTSliceToVolumeRegistration.h"

EBBRT_PUBLISH_TYPE(, EbbRTSliceToVolumeRegistration);

void EbbRTSliceToVolumeRegistration::Print(ebbrt::Messenger::NetworkId nid,
                                           const char* str) {
  auto len = strlen(str) + 1;
  auto buf = ebbrt::MakeUniqueIOBuf(len);
  snprintf(reinterpret_cast<char*>(buf->MutData()), len, "%s", str);
  SendMessage(nid, std::move(buf));
}

void EbbRTSliceToVolumeRegistration::ReceiveMessage(
    ebbrt::Messenger::NetworkId nid, std::unique_ptr<ebbrt::IOBuf>&& buffer) {
  /***********************************
   * std::string(reinterpret_cast<const char*> buffer->Data(),
   *buffer->Length())
   *
   * Using the above code, I had to first do strcmp(output.c_str(), ...) to
   *ensure
   * it matched the input string.
   * Direct comparison using "==" seems to be working when I don't pass in the
   *Length() as second arg
   **********************************/
  auto output = std::string(reinterpret_cast<const char*>(buffer->Data()));
  std::cout << "Received msg: " << nid.ToString() << ": " << output << "\n";

  // stores the received nid as a string of bytes
  if (bmnid.empty()) {
    bmnid = nid.ToBytes();
    std::cout << "isEmpty()? " << bmnid.empty() << std::endl;

    //*****
    // this SetValue() sets up the promise so that the waitReceive()
    // future is fulfilled and starts the next computation in AppMain()
    // this call ensures that we don't call SendMessage until we know
    // the bm node has been initialized
    mypromise.SetValue();
  }

  if (output[0] == 'A') {
    std::cout << "Parsing it back" << std::endl;

    std::istringstream ifs;
    ifs.str(output.substr(2, output.length() - 2));
    boost::archive::text_iarchive ia(ifs);

    ia& mym;

    // reactivate context in runJob2()
    ebbrt::event_manager->ActivateContext(std::move(*emec));
  } 
  else if (output[0] == 'B') 
  {
      std::cout << "Parsing it back" << std::endl;

      std::istringstream ifs;
      ifs.str(output.substr(2, output.length() - 2));
      boost::archive::text_iarchive ia(ifs);

      ia & t & trans & regCert; 
            
      // reactivate context in runJob2()
      ebbrt::event_manager->ActivateContext(std::move(*emec));
    }
}

void EbbRTSliceToVolumeRegistration::runJob() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  irtkMatrix m1(2, 2);
  irtkMatrix m2(2, 2);
  // irtkMatrix m3;

  initMatrix(m1, 2.0, 2, 2);
  initMatrix(m2, 3.0, 2, 2);

  std::ostringstream ofs;
  boost::archive::text_oarchive oa(ofs);
  oa& m1& m2;

  std::string ts = "A " + ofs.str();

  // sending serialized vector to bm
  ebbrt::Messenger::NetworkId bmnidd;
  Print(bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
                         bmnid.size()),
        ts.c_str());
  std::cout << "Sent vector " << std::endl;

  emec = &context;
  ebbrt::event_manager->SaveContext(*emec);

  // this should execute when the saved context is loaded again
  std::cout << "Received matrix: " << std::endl;
  mym.Print();
  ebbrt::active_context->io_service_.stop();
}

void EbbRTSliceToVolumeRegistration::slicetovolume() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  irtkImageAttributes attr = reconstructor->_reconstructed.GetImageAttributes();
  size_t max_slices = reconstructor->_slices.size();
  size_t inputIndex = 0;
  
  for (inputIndex = 0; inputIndex != max_slices; inputIndex++) 
  {
      // serialization
      ebbrt::Messenger::NetworkId bmnidd;
      std::ostringstream ofs;
      boost::archive::text_oarchive oa(ofs);
      
      oa& reconstructor->_slices[inputIndex] & attr & reconstructor->_reconstructed
	  & reconstructor->_transformations[inputIndex];
      
      std::string ts = "B " + ofs.str();
      Print(
	  bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
			   bmnid.size()),
	  ts.c_str());

      emec = &context;
      ebbrt::event_manager->SaveContext(*emec);

      reconstructor->_slices[inputIndex] = t;
      reconstructor->_transformations[inputIndex] = trans;
      reconstructor->_slices_regCertainty[inputIndex] = regCert;

      std::cout << "EbbRTSliceToVolumeRegistration done " << std::endl;
      /*
    irtkImageRigidRegistrationWithPadding registration;
    irtkGreyPixel smin, smax;
    irtkGreyImage target;
    irtkRealImage slice, w, b, t;
    irtkResamplingWithPadding<irtkRealPixel> resampling(attr._dx, attr._dx,
                                                        attr._dx, -1);

    // t = reconstructor->_slices[inputIndex];

    // serialization
    ebbrt::Messenger::NetworkId bmnidd;
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);

    oa& reconstructor->_slices[inputIndex] & attr;

    std::string ts = "B " + ofs.str();
    Print(
        bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
                         bmnid.size()),
        ts.c_str());
    std::cout << "Sent vector " << std::endl;
    emec = &context;
    ebbrt::event_manager->SaveContext(*emec);
    
    
    // resampling.SetInput(&reconstructor->_slices[inputIndex]);
    // resampling.SetOutput(&t);
    // resampling.Run();
    
    //reconstructor->_slices[inputIndex] = t;
    target = t;

    target.GetMinMax(&smin, &smax);

    if (smax > -1) {
      // put origin to zero
      irtkRigidTransformation offset;
      irtkReconstruction::ResetOrigin(target, offset);

      irtkMatrix mo = offset.GetMatrix();
      irtkMatrix m = reconstructor->_transformations[inputIndex].GetMatrix();
      m = m * mo;
      reconstructor->_transformations[inputIndex].PutMatrix(m);

      irtkGreyImage source = reconstructor->_reconstructed;
      registration.SetInput(&target, &source);
      registration.SetOutput(&reconstructor->_transformations[inputIndex]);
      registration.GuessParameterSliceToVolume();
      registration.SetTargetPadding(-1);
      registration.Run();

      reconstructor->_slices_regCertainty[inputIndex] =
          registration.last_similarity;

      // undo the offset
      mo.Invert();
      m = reconstructor->_transformations[inputIndex].GetMatrix();
      m = m * mo;
      reconstructor->_transformations[inputIndex].PutMatrix(m);
    }
    }*/
  }
}
