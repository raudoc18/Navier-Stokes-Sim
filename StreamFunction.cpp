//
// Created by Dominik on 25.09.26.
//

#include "Streamfunction.h";
#include "constants.h";

typedef datastruct::Matrix<double> Matrix;

void StreamFunction::computeStreamFunc() {
    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize + 1;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

        for (int j = start_idx; j < end_idx; j++) {
            for (int i = 1; i < nx + 1; i++) {
                streamFunc(j, i) = streamFunc(j, i - 1) + u(j, i) * dx;
            }
        }
    });
    printf("foo");
}

void StreamFunction::extractLines() {
    double max = streamFunc.data()[streamFunc.maxIdx()];
    double min = streamFunc.data()[streamFunc.minIdx()];


}
