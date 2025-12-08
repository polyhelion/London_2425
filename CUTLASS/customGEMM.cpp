
// GEMM kernel with custom custom tiling config 

using GemmKernel = cutlass::gemm::kernel::DefaultGemm<
    float, cutlass::layout::RowMajor,   // A
    float, cutlass::layout::RowMajor,   // B
    float, cutlass::layout::RowMajor,   // C
    float,                              // Accumulator
    cutlass::arch::OpClassTensorOp,     // Use tensor cores
    cutlass::arch::Sm80,                // GPU architecture
    cutlass::gemm::GemmShape<128, 128, 32>,   // Threadblock tile
    cutlass::gemm::GemmShape<64, 64, 32>,     // Warp tile
    cutlass::gemm::GemmShape<16, 8, 8>,       // Instruction (tensor core) tile
    cutlass::epilogue::thread::LinearCombination<
        float, 1, float, float>,        // Epilogue operation (C = alpha * A*B + beta*C)
    cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>, // Mapping
    3                                  // Stages in pipeline (for double-buffering)
>::GemmKernel;

// device-level operator
using Gemm = cutlass::gemm::device::GemmUniversalAdapter<GemmKernel>;


