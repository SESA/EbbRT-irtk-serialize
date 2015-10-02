//          Copyright Boston University SESA Group 2013 - 2014.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)
#include "Printer.h"

#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include <string.h>
#include <vector>

#include <boost/serialization/vector.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <ebbrt/GlobalIdMap.h>
#include <ebbrt/UniqueIOBuf.h>

EBBRT_PUBLISH_TYPE(, Printer);

Printer::Printer() : ebbrt::Messagable<Printer>(kPrinterEbbId) {}

ebbrt::Future<void> Printer::Init() {
  return ebbrt::global_id_map->Set(
      kPrinterEbbId, ebbrt::messenger->LocalNetworkId().ToBytes());
}

void Printer::Print(ebbrt::Messenger::NetworkId nid, const char* str) {
  auto len = strlen(str) + 1;
  auto buf = ebbrt::MakeUniqueIOBuf(len);
  snprintf(reinterpret_cast<char*>(buf->MutData()), len, "%s", str);
  SendMessage(nid, std::move(buf));
}

void Printer::ReceiveMessage(ebbrt::Messenger::NetworkId nid,
                             std::unique_ptr<ebbrt::IOBuf>&& buffer) {
  /***********************************
   * std::string(reinterpret_cast<const char*> buffer->Data(), buffer->Length())
   *
   * Using the above code, I had to first do strcmp(output.c_str(), ...) to
   *ensure
   * it matched the input string.
   * Direct comparison using "==" seems to be working when I don't pass in the
   *Length() as second arg
   **********************************/
  auto output = std::string(reinterpret_cast<const char*>(buffer->Data()));
  std::this_thread::sleep_for(std::chrono::seconds(1));

  // serialize a vector
  std::vector<int> v;
  v.push_back(1);
  v.push_back(2);
  v.push_back(3);
  std::ostringstream ofs;
  boost::archive::text_oarchive oa(ofs);
  oa& v;

  if (output == "ping") {
    std::cout << nid.ToString() << ": " << output << "\n";
    Print(nid, ofs.str().c_str());
  } else if (strcmp(output.c_str(), "pong") == 0) {
    std::cout << nid.ToString() << ": " << output << "\n";
    Print(nid, "ping");
  }
  
  std::this_thread::sleep_for(std::chrono::seconds(1));
  std::cout << "Shutting down ..." << std::endl;
  
  // referencing active_context defined in EbbRT/hosted/src/include/ebbrt/Context.h
  // gets this current ebbs active context
  ebbrt::active_context->io_service_.stop();
}
