#include "EbbRTCoeffInit.h"

EBBRT_PUBLISH_TYPE(, EbbRTCoeffInit);

struct membuf : std::streambuf {
  membuf(char* begin, char* end) { this->setg(begin, begin, end); }
};

void EbbRTCoeffInit::Print(ebbrt::Messenger::NetworkId nid, const char* str) {
  auto len = strlen(str) + 1;
  auto buf = ebbrt::MakeUniqueIOBuf(len);
  snprintf(reinterpret_cast<char*>(buf->MutData()), len, "%s", str);

  std::cout << "runJob2(): length of sent iobuf: "
            << buf->ComputeChainDataLength() << " bytes" << std::endl;

  SendMessage(nid, std::move(buf));
}

void EbbRTCoeffInit::ReceiveMessage(ebbrt::Messenger::NetworkId nid,
                                    std::unique_ptr<ebbrt::IOBuf>&& buffer) {
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
  std::cout << "Received ip: " << nid.ToString() << std::endl;

  //  std::cout << "size: " << vNids.size() << std::endl;
  // stores the received nid as a string of bytes
  // if (bmnid.empty()) {
  // if ((int)vNids.size() != numNodes) {
  //      std::cout << "size: " << vNids.size() << std::endl;
  // ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
  // char *t = (char *)(dp.Get(buffer->ComputeChainDataLength()));
  // std::string st(t);
  // std::cout << "chain length: " << buffer->ComputeChainDataLength() << "\nst:
  // " << st << std::endl;

  // bmnid = nid.ToBytes();
  //  std::cout << "Received msg: " << nid.ToString() << ": " << output << "\n";
  //  vNids.push_back(nid.ToBytes());

  // this SetValue() sets up the promise so that the waitReceive()
  // future is fulfilled and starts the next computation in AppMain()
  // this call ensures that we don't call SendMessage until we know
  // the bm node has been initialized
  //  if ((int)vNids.size() == numNodes) {
  //	  mypromise.SetValue();
  //    }
  //}

  if (output[0] == 'A') {
    std::cout << "Parsing it back" << std::endl;

    std::istringstream ifs;
    ifs.str(output.substr(2, output.length() - 2));
    boost::archive::text_iarchive ia(ifs);

    ia& mym;

    // reactivate context in runJob2()
    ebbrt::event_manager->ActivateContext(std::move(*emec));
  } else if (output[0] == 'C') {
    std::cout << "Parsing it back" << std::endl;

    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};

    boost::archive::text_iarchive ia(stream);

    //      std::istringstream ifs;
    //    ifs.str(output.substr(2, output.length() - 2));
    // boost::archive::text_iarchive ia(ifs);

    ia& reconstructor->_slices& reconstructor->_transformations& reconstructor
        ->_volcoeffs& reconstructor->_slice_inside_cpu;

    // ia & _slice & _transformations & slicecoeffs & slice_inside;

    // reactivate context in runJob2()
    ebbrt::event_manager->ActivateContext(std::move(*emec));
  } else if (output[0] == 'D') {
    std::cout << "FE received D" << std::endl;
    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};

    boost::archive::text_iarchive ia(stream);
    ia& _slice;

    std::cout << "Parsing it back, received: "
              << buffer->ComputeChainDataLength() << " bytes" << std::endl;
    recv_counter++;

    if (recv_counter == numNodes) {
      // reactivate context in runJob2()
      ebbrt::event_manager->ActivateContext(std::move(*emec));
    }
  } else if (output[0] == 'G') {
    std::cout << "FE received D" << std::endl;
    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};
    std::vector<int> rvec;

    boost::archive::text_iarchive ia(stream);
    ia& rvec;

    std::cout << "Parsing it back, received: "
              << buffer->ComputeChainDataLength() << " bytes" << std::endl;
    std::cout << "rvec.size() = " << rvec.size() << std::endl;
    
    recv_counter++;
    std::cout << "recv_counter++" << recv_counter << std::endl;

    if (recv_counter == numNodes) {
	std::cout << "recv_counter == numNodes" << recv_counter << std::endl;
      // reactivate context
      ebbrt::event_manager->ActivateContext(std::move(*emec));
    }
  } else if (output[0] == 'E') {
    ebbrt::IOBuf::DataPointer dp = buffer->GetDataPointer();
    char* t = (char*)(dp.Get(buffer->ComputeChainDataLength()));
    membuf sb{t + 2, t + buffer->ComputeChainDataLength()};
    std::istream stream{&sb};
    boost::archive::text_iarchive ia(stream);

    std::cout << "Parsing it back, received: "
              << buffer->ComputeChainDataLength() << " bytes" << std::endl;

    int start, end;
    ia& start& end;

    std::cout << "start: " << start << " end: " << end << std::endl;

    /*    for(int k = start; k < end; k++) { ia & reconstructor->_slices[k]; }
        std::cout << "deserialize _slices" << std::endl;
        for(int k = start; k < end; k++) { ia &
       reconstructor->_transformations[k]; }
        std::cout << "deserialize _transformations" << std::endl;
        for(int k = start; k < end; k++) { ia & reconstructor->_volcoeffs[k]; }
        std::cout << "deserialize _volcoeffs" << std::endl;*/

    std::vector<bool> _slice_inside_cpu_sub;
    std::vector<irtkRealImage> _slices_sub;
    std::vector<irtkRigidTransformation> _transformations_sub;
    std::vector<SLICECOEFFS> _volcoeffs_sub;

    ia& _slices_sub& _transformations_sub& _volcoeffs_sub&
    _slice_inside_cpu_sub;
    std::cout << "deserialize _slices " << _slices_sub.size() << std::endl;
    std::cout << "deserialize _transformations" << _transformations_sub.size() << std::endl;
    std::cout << "deserialize _volcoeffs" << _volcoeffs_sub.size() << std::endl;

    int c = 0;
    
    for (int k = start; k < end; k++) {
	reconstructor->_slice_inside_cpu[k] = _slice_inside_cpu_sub[c];
	reconstructor->_slices[k] = _slices_sub[c];
	reconstructor->_transformations[k] = _transformations_sub[c];
	reconstructor->_volcoeffs[k] = _volcoeffs_sub[c];
	c++;
    }

    std::cout << "deserialize _slice_inside_cpu" << std::endl;

    //    ia& reconstructor->_slices& reconstructor->_transformations&
    // reconstructor
    //      ->_volcoeffs& reconstructor->_slice_inside_cpu;

    recv_counter++;
    
    if (recv_counter == numNodes) {
      // reactivate context
      ebbrt::event_manager->ActivateContext(std::move(*emec));
    }
  } else {
    std::cout << "Received msg: " << nid.ToString() << ": " << output << "\n";
  }
}

void EbbRTCoeffInit::runJob3() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  // serialization
  std::ostringstream ofs;
  boost::archive::text_oarchive oa(ofs);
  oa& reconstructor->_slices[0];

  ebbrt::Messenger::NetworkId bmnidd;
  std::string ts = "E " + ofs.str();

  Print(bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
                         bmnid.size()),
        ts.c_str());

  emec = &context;
  ebbrt::event_manager->SaveContext(*emec);

  std::cout << "Received back " << std::endl;
  reconstructor->_slices[0] = _slice;

  bmnid.clear();
  ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::runJob(int size) {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  std::vector<int> vec; //{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  for(int i = 0; i < size; i ++)
  {
      vec.push_back(i);
  }

  for (int i = 0; i < (int)nids.size(); i++) {
    int start, end, factor;
    factor = (int)ceil(vec.size() / (float)numNodes);
    start = i * factor;
    end = i * factor + factor;
    end = ((size_t)end > vec.size()) ? vec.size() : end;

    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);

    oa& start& end;

    for (int j = start; j < end; j++) {
      oa& vec[j];
    }

    std::string ts = "G " + ofs.str();

    std::cout << "Sending to .. " << nids[i].ToString() << " " << start << " "
              << end << std::endl;
    Print(nids[i], ts.c_str());
  }

  emec = &context;
  ebbrt::event_manager->SaveContext(*emec);

  std::cout << "Received back " << std::endl;
  ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::runJob2() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  // serialization
  ebbrt::Messenger::NetworkId bmnidd;
  // std::string ts = "ping";
  std::string ts = "D " + reconstructedstr;

  // Print(nid, ts.c_str());
  for (int i = 0; i < (int)nids.size(); i++) {
    std::cout << "Sending to .. " << nids[i].ToString() << std::endl;
    //	Print(bmnidd.FromBytes(reinterpret_cast<const unsigned
    // char*>(vNids[i].c_str()),
    //			       vNids[i].size()), ts.c_str());
    Print(nids[i], ts.c_str());
  }

  std::cout << "Sending data" << std::endl;
  emec = &context;
  ebbrt::event_manager->SaveContext(*emec);

  std::cout << "Received back " << std::endl;
  ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::coeffinit() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  irtkImageAttributes attr = reconstructor->_reconstructed.GetImageAttributes();
  size_t max_slices = reconstructor->_slices.size();
  size_t inputIndex = 0;

  std::cout << "Inside coeffinit() " << max_slices << " " << inputIndex
            << std::endl;

  for (int i = 0; i < (int)nids.size(); i++) {

    // serialization
    // ebbrt::Messenger::NetworkId bmnidd;
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);

    std::cout << "Before serialize" << std::endl;

    oa& reconstructor->_reconstructed& reconstructor
        ->_quality_factor& reconstructor->_mask& max_slices& reconstructor
        ->_slices& reconstructor->_transformations& reconstructor
        ->_volcoeffs& reconstructor->_slice_inside_cpu;

    /*for (inputIndex = 0; inputIndex != max_slices; inputIndex++)
    {
        oa & reconstructor->_slices[inputIndex]
            & reconstructor->_transformations[inputIndex];
            }*/

    std::cout << "after serialize" << std::endl;

    std::string ts = "C " + ofs.str();
    // Print(bmnidd.FromBytes(reinterpret_cast<const unsigned
    // char*>(bmnid.c_str()),
    //                     bmnid.size()),
    //    ts.c_str());
    
    std::cout << "Sending to .. " << nids[i].ToString()
              << " size: " << ts.length() << std::endl;
    
    Print(nids[i], ts.c_str());
  }

  std::cout << "Sending data" << std::endl;
  emec = &context;
  ebbrt::event_manager->SaveContext(*emec);

  std::cout << "Received back " << std::endl;
  // reconstructor->_slices[inputIndex] = _slice;
  // reconstructor->_transformations[inputIndex] = _transformations;
  // reconstructor->_volcoeffs[inputIndex] = slicecoeffs;
  // reconstructor->_slice_inside_cpu[inputIndex] = slice_inside;

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
          bmnidd.FromBytes(reinterpret_cast<const unsigned
  char*>(bmnid.c_str()),
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

  std::cout << "Inside coeffinit() " << max_slices << " " << inputIndex
            << std::endl;

  // serialization
  ebbrt::Messenger::NetworkId bmnidd;
  std::ostringstream ofs;
  boost::archive::text_oarchive oa(ofs);

  std::cout << "Before serialize" << std::endl;

  oa& reconstructor->_reconstructed& reconstructor
      ->_quality_factor& reconstructor->_mask& max_slices& reconstructor
      ->_slices& reconstructor->_transformations& reconstructor
      ->_volcoeffs& reconstructor->_slice_inside_cpu;

  std::cout << "after serialize" << std::endl;

  std::string ts = "F " + ofs.str();

  for (int i = 0; i < (int)nids.size(); i++) {
    std::cout << "Sending to .. " << nids[i].ToString() << std::endl;
    Print(nids[i], ts.c_str());
  }
  //    Print(
  //	bmnidd.FromBytes(reinterpret_cast<const unsigned char*>(bmnid.c_str()),
  // bmnid.size()),
  //	ts.c_str()
  //	);

  std::cout << "Sending data" << std::endl;
  emec = &context;
  ebbrt::event_manager->SaveContext(*emec);

  std::cout << "Received back " << std::endl;
  std::cout << "EbbRTCoeffInit done " << std::endl;
  ebbrt::active_context->io_service_.stop();
}

void EbbRTCoeffInit::coeffinitParallel2() {
  // get the event manager context and save it
  ebbrt::EventManager::EventContext context;

  irtkImageAttributes attr = reconstructor->_reconstructed.GetImageAttributes();
  size_t max_slices = reconstructor->_slices.size();
  size_t inputIndex = 0;

  std::cout << "Inside coeffinitParallel2() " << max_slices << " " << inputIndex
            << std::endl;

  int start, end, factor;

  // all sizes are the same
  std::cout << "_slices.size() " << reconstructor->_slices.size() << std::endl;
  std::cout << "_transformations.size() "
            << reconstructor->_transformations.size() << std::endl;
  std::cout << "__volcoeffs.size() " << reconstructor->_volcoeffs.size()
            << std::endl;
  std::cout << "___slice_inside_cpu.size() "
            << reconstructor->_slice_inside_cpu.size() << std::endl;

  for (int i = 0; i < (int)nids.size(); i++) {
    // serialization
    std::ostringstream ofs;
    boost::archive::text_oarchive oa(ofs);

    std::cout << "Before serialize" << std::endl;

    //_slices
    factor = (int)ceil(max_slices / (float)numNodes);
    start = i * factor;
    end = i * factor + factor;
    end = ((size_t)end > max_slices) ? max_slices : end;

    oa& start& end;

    /*    for (int j = start; j < end; j++) {
          oa& reconstructor->_slices[j];
        }

        for (int j = start; j < end; j++) {
          oa& reconstructor->_transformations[j];
        }

        for (int j = start; j < end; j++) {
          oa& reconstructor->_volcoeffs[j];
          }*/

    std::vector<bool> _slice_inside_cpu_sub;
    std::vector<irtkRealImage> _slices_sub;
    std::vector<irtkRigidTransformation> _transformations_sub;
    std::vector<SLICECOEFFS> _volcoeffs_sub;

    for (int j = start; j < end; j++) {
      _slice_inside_cpu_sub.push_back(reconstructor->_slice_inside_cpu[j]);
      _slices_sub.push_back(reconstructor->_slices[j]);
      _volcoeffs_sub.push_back(reconstructor->_volcoeffs[j]);
      _transformations_sub.push_back(reconstructor->_transformations[j]);
    }

    oa& _slices_sub& _transformations_sub& _volcoeffs_sub&
    _slice_inside_cpu_sub;
    // oa& reconstructor->_slice_inside_cpu;

    oa& reconstructor->_reconstructed& reconstructor
        ->_quality_factor& reconstructor->_mask& max_slices;

    std::string ts = "E " + ofs.str();

    std::cout << "Sending to .. " << nids[i].ToString()
              << " size: " << ts.length() << std::endl;
    Print(nids[i], ts.c_str());
  }

  std::cout << "Saving context " << std::endl;

  emec = &context;
  ebbrt::event_manager->SaveContext(*emec);

  std::cout << "Received back " << std::endl;
  std::cout << "EbbRTCoeffInit done " << std::endl;
  ebbrt::active_context->io_service_.stop();
}
