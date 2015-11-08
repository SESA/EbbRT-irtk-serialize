#include "EbbRTCoeffInit.h"

EBBRT_PUBLISH_TYPE(, EbbRTCoeffInit);
//EBBRT_PUBLISH_TYPE(, Printer);

struct membuf : std::streambuf {
    membuf(char* begin, char* end) {
	this->setg(begin, begin, end);
    }
};

void EbbRTCoeffInit::Print(ebbrt::Messenger::NetworkId nid,
                                           const char* str) {
  auto len = strlen(str) + 1;
  auto buf = ebbrt::MakeUniqueIOBuf(len);
  snprintf(reinterpret_cast<char*>(buf->MutData()), len, "%s", str);

  std::cout << "runJob2(): length of sent iobuf: " << buf->ComputeChainDataLength() << " bytes" << std::endl;

  SendMessage(nid, std::move(buf));
}

void EbbRTCoeffInit::ReceiveMessage(
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
  //std::cout << "Received msg: " << nid.ToString() << ": " << output << "\n";

  // stores the received nid as a string of bytes
  if (bmnid.empty()) {
  //if ((int)vNids.size() != numNodes) {
      //std::cout << "size: " << vNids.size() << std::endl;
      //ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
      //char *t = (char *)(dp.Get(buffer->ComputeChainDataLength()));
      //std::string st(t);
      //std::cout << "chain length: " << buffer->ComputeChainDataLength() << "\nst: " << st << std::endl;
      
      bmnid = nid.ToBytes();
      //vNids.push_back(nid.ToBytes());
      
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
  else if (output[0] == 'C') 
  {
      std::cout << "Parsing it back" << std::endl;

      ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
      char *t = (char *)(dp.Get(buffer->ComputeChainDataLength()));
      membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
      std::istream stream{&sb};
            
      boost::archive::text_iarchive ia(stream);

//      std::istringstream ifs;
      //    ifs.str(output.substr(2, output.length() - 2));
      //boost::archive::text_iarchive ia(ifs);

      ia & reconstructor->_slices
	  & reconstructor->_transformations
	  & reconstructor->_volcoeffs
	  & reconstructor-> _slice_inside_cpu;
      
      //ia & _slice & _transformations & slicecoeffs & slice_inside;
      
      // reactivate context in runJob2()
      ebbrt::event_manager->ActivateContext(std::move(*emec));
  }
  else if (output[0] == 'D')
  {
      std::cout << "FE received D" << std::endl;
      ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
      char *t = (char *)(dp.Get(buffer->ComputeChainDataLength()));
      membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
      std::istream stream{&sb};
            
      boost::archive::text_iarchive ia(stream);
      ia & _slice;

      std::cout << "Parsing it back, received: " << buffer->ComputeChainDataLength() << " bytes" << std::endl;

      // reactivate context in runJob2()
      ebbrt::event_manager->ActivateContext(std::move(*emec));
  }
  else if (output[0] == 'E')
  {
      ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
      char *t = (char *)(dp.Get(buffer->ComputeChainDataLength()));
      membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
      std::istream stream{&sb};
            
      boost::archive::text_iarchive ia(stream);
      ia & _slice;

      std::cout << "Parsing it back, received: " << buffer->ComputeChainDataLength() << " bytes" << std::endl;

      //ia & _slice & _transformations & slicecoeffs & slice_inside;
      
      // reactivate context in runJob2()
      ebbrt::event_manager->ActivateContext(std::move(*emec));
  }
}

void EbbRTCoeffInit::runJob3()
{
    // get the event manager context and save it
    ebbrt::EventManager::EventContext context;

    // serialization
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
    oa & reconstructor->_slices[0];
    
    ebbrt::Messenger::NetworkId bmnidd;
    std::string ts = "E " + ofs.str();
    
    Print(bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
			   bmnid.size()), ts.c_str());
    
    emec = &context;
    ebbrt::event_manager->SaveContext(*emec);
    
    std::cout << "Received back " <<std::endl;
    reconstructor->_slices[0] = _slice;

    bmnid.clear();
    ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::runJob() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;
  int i = 0;
  
  for(i=0;i<2;i++)
  {
      irtkMatrix m1(2, 2);
      irtkMatrix m2(2, 2);
      // irtkMatrix m3;

      initMatrix(m1, 2.0, 2, 2);
      initMatrix(m2, 3.0, 2, 2);
      irtkRealImage test;
      
      std::cout << "before" <<std::endl;
      
      std::ostringstream ofs;
      boost::archive::text_oarchive oa(ofs);
      oa & m1 & m2 & test;
      
      std::cout << "after" <<std::endl;

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
  }

  ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::runJob2()
{    
    // get the event manager context and save it
    ebbrt::EventManager::EventContext context;

    // serialization
    ebbrt::Messenger::NetworkId bmnidd;
    //std::string ts = "D abc";
    std::string ts = "D " + reconstructedstr;
  
    Print(bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
			   bmnid.size()), ts.c_str());

    std::cout << "Sending data" << std::endl;
    emec = &context;
    ebbrt::event_manager->SaveContext(*emec);
    
    std::cout << "Received back " <<std::endl;
    ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::coeffinit() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  irtkImageAttributes attr = reconstructor->_reconstructed.GetImageAttributes();
  size_t max_slices = reconstructor->_slices.size();
  size_t inputIndex = 0;
  
  std::cout << "Inside coeffinit() " << max_slices 
	    << " " << inputIndex << std::endl;

  // serialization
  ebbrt::Messenger::NetworkId bmnidd;
  std::ostringstream ofs;
  boost::archive::text_oarchive oa(ofs);
  
  std::cout << "Before serialize" << std::endl;

  oa & reconstructor->_reconstructed
      & reconstructor->_quality_factor 
      & reconstructor->_mask
      & max_slices
      & reconstructor->_slices
      & reconstructor->_transformations
      & reconstructor->_volcoeffs
      & reconstructor->_slice_inside_cpu;
  
  /*for (inputIndex = 0; inputIndex != max_slices; inputIndex++) 
  {
      oa & reconstructor->_slices[inputIndex]
	  & reconstructor->_transformations[inputIndex];
	  }*/
  
  std::cout << "after serialize" << std::endl;
      
  std::string ts = "C " + ofs.str();
  Print(
      bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()), bmnid.size()),
      ts.c_str()
  );
  
  std::cout << "Sending data" << std::endl;
  emec = &context;
  ebbrt::event_manager->SaveContext(*emec);
  
  std::cout << "Received back " <<std::endl;
  //reconstructor->_slices[inputIndex] = _slice;
  //reconstructor->_transformations[inputIndex] = _transformations;
  //reconstructor->_volcoeffs[inputIndex] = slicecoeffs;
  //reconstructor->_slice_inside_cpu[inputIndex] = slice_inside;
  
  std::cout << "EbbRTCoeffInit done " << std::endl;

  /*for (inputIndex = 0; inputIndex != max_slices; inputIndex++) 
  {
      // serialization
      ebbrt::Messenger::NetworkId bmnidd;
      std::ostringstream ofs;
      boost::archive::text_oarchive oa(ofs);

      std::cout << "Before serialize" << std::endl;
      
      oa & reconstructor->_slices[inputIndex]
	  & reconstructor->_reconstructed
	  & reconstructor->_quality_factor 
	  & reconstructor->_transformations[inputIndex]
	  & reconstructor->_mask;
      
      std::cout << "after serialize" << std::endl;
      
      std::string ts = "C " + ofs.str();
      Print(
	  bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
			   bmnid.size()),
	  ts.c_str());

      std::cout << "Sending data" << std::endl;
      emec = &context;
      ebbrt::event_manager->SaveContext(*emec);

      std::cout << "Received back " <<std::endl;
      reconstructor->_slices[inputIndex] = _slice;
      reconstructor->_transformations[inputIndex] = _transformations;
      reconstructor->_volcoeffs[inputIndex] = slicecoeffs;
      reconstructor->_slice_inside_cpu[inputIndex] = slice_inside;
      
      std::cout << "EbbRTCoeffInit done " << std::endl;
  }
  */
  ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::coeffinitParallel() {
    // get the event manager context and save it
    ebbrt::EventManager::EventContext context;

    irtkImageAttributes attr = reconstructor->_reconstructed.GetImageAttributes();
    size_t max_slices = reconstructor->_slices.size();
    size_t inputIndex = 0;
  
    std::cout << "Inside coeffinit() " << max_slices 
	      << " " << inputIndex << std::endl;

    // serialization
    ebbrt::Messenger::NetworkId bmnidd;
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);
  
    std::cout << "Before serialize" << std::endl;

    oa & reconstructor->_reconstructed
	& reconstructor->_quality_factor 
	& reconstructor->_mask
	& max_slices
	& reconstructor->_slices
	& reconstructor->_transformations
	& reconstructor->_volcoeffs
	& reconstructor->_slice_inside_cpu;

    std::cout << "after serialize" << std::endl;
      
    std::string ts = "F " + ofs.str();
    Print(
	bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()), bmnid.size()),
	ts.c_str()
	);
  
    std::cout << "Sending data" << std::endl;
    emec = &context;
    ebbrt::event_manager->SaveContext(*emec);
  
    std::cout << "Received back " <<std::endl;
    std::cout << "EbbRTCoeffInit done " << std::endl;
    ebbrt::active_context->io_service_.stop();
}
