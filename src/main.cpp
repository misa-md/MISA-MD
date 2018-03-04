//
// Created by gensh(genshenchu@gmail.com)  on 2017/4/15.
//

#include "crystal_md.h"

int main(int argc, char **argv) {
    // app's lifecycle here.
    auto *app = new crystalMD(argc, argv);
    if ((app->initialize())) {
        if (app->prepare()) {
            app->run();
            app->destroy();
        }
    }
    app->detach();
    return 0;
}
