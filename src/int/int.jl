export int, DIRECT, EXCHANGE, PAIR

@enum InterAction DIRECT EXCHANGE PAIR

include("./_report_rhok.jl")

# preprocessing is the key step.
include("./preprocessing/preprocessing.jl")
include("./postprocessing/postprocessing.jl")

include("./direct/direct.jl")
include("./exchange/exchange.jl")
include("./pair/pair.jl")

include("./woeps/woeps.jl")
include("./wieps/wieps.jl")

# 关于BLAS多线程：默认的OpenBLAS的行为是其维护自己的一个线程池，对外表现串行(但是我测试的结果表明Julia多线程内部的BLAS操作会独立进行，并且应该是单线程，并不会排队等待OpenBLAS线程池)；
# MKL则是会为每个Julia线程分发n个线程，最终有可能导致线程超发。
# 比较通用且合适的方式为：避免使用双重多线程；尽量使用Julia多线程，并明确关闭BLAS多线程；个别大型计算任务使用BLAS多线程一次性处理。
# 关于OpenBLAS的多线程，具体是否启用依赖于计算任务的大小，其会自动判断是否能够抵消线程调用的开销。
# 关于FFTW，如果单个FFTW任务过大(2D及以下)，内存带宽成为主要瓶颈，盲目增加线程数可能并不能显著降低计算总时长。
# FFTW预计会成为主要瓶颈，优化路径为手动进行多次1DFFT来替换完整的3DFFT，如此可以在Julia层级进行更精细的控制，并避免内存带宽的限制。
# 预计需要使用permutedims!来避免跨步访问，保证FFT进行的始终是连续内存。

# 一个案例，72核，ρR的计算过程中大概只占用了10核左右。
# 需要将ρk和ρR的计算结合起来，将任务粒度划分更细。
"""
    int(ρk::AbstractRhok; kwargs...)
    int(epsilon::DielectricPWBasis, ρk::AbstractRhok; kwargs...)

Evaluate coulomb interaction in Wannier basis, with or without screening.
Only support 3D and 2D system at present.

- `pair`::Tuple{<:Integer,<:Integer,<:AbstractVector{<:Integer}}
    A tuple `(i, j, R)` representing a pair formed by the Wannier basis at `i0` and the Wannier basis at `jR`.
- `pair_uplo`::Symbol = :all
    Whether to reduce the computational cost by exploiting the Hermiticity of the inter-pair interactions. Supported options are `:all`, `:U`, and `:L`.
- `pair_atol`::Real = sqrt(eps(Float64))
    Truncation threshold for the final computed inter-pair interactions. Values below this threshold are discarded and are not saved.
- `coulomb`
    Truncated coulomb potential in plane-wave basis, no need to set up for 3D and 2D system, and have to set up for 1D and 0D system.
- `outfolder`::AbstractString = "./"
    Directory for saving results in JLD2 format. Files are suffixed with either `bare` or an integer, representing results at different frequency points and matching the frequency points of the dielectric function one-to-one.

A log file named `WannierInterAction.log` is written to `outfolder`; 
it can be inspected during the run to track which tasks have completed and to estimate the time remaining.

See also [`Rhok`](@ref), [`PlaneWaveCoulomb`](@ref).
"""
function int(ρk::AbstractRhok;
    ρRnchunk::Integer=Threads.nthreads(), ρRbuffer::Symbol=:auto,
    iqidx=:default, fileformat::Symbol=:jld2,
    outfolder::AbstractString="./", logfile::AbstractString=joinpath(outfolder, "./WannierInterAction.log"),
    interaction::AbstractVector{InterAction}=[DIRECT, EXCHANGE, PAIR],
    coulomb::Union{PlaneWaveCoulomb,Nothing}=nothing,
    multipolecorrection::Union{AbstractMultipoleCorrection,Nothing}=nothing,
    wingaverage::Bool=true,
    grid_atol::Real=1e-6, exchange_range=(3, 3, 3),
    pair=:default, pair_uplo=:all, pair_range=:default, pair_atol::Real=sqrt(eps(Float64)),
)

    # logio
    mkpath(dirname(logfile))
    logio = get_logio(logfile)
    original_stdout = stdout
    redirect_stdout(logio)

    # used for debugging
    if iqidx == :default
        iqidx = collect(1:ρk.nq)
    else
        iqidx = Int.(collect(iqidx))
        if any(i -> i < 1 || i > ρk.nq, iqidx)
            throw(ArgumentError("iqidx must be between 1 and $(ρk.nq)"))
        end
    end

    interaction = Set(interaction)
    @assert length(interaction) > 0 "Empty `interaction`."

    # check some input parameters that cannot be automatically pre-checked by the program.
    @assert fileformat ∈ [:jld2, :txt] "Wrong `fileformat`!"

    interaction_para = (
        interaction=interaction,
        wingaverage=wingaverage,
        grid_atol=grid_atol,
        exchange_range=exchange_range,
        pair=pair,
        pair_uplo=pair_uplo,
        pair_range=pair_range,
        pair_atol=pair_atol,
        ρR_buffer_level=ρRbuffer,
        JuliaTask_nchunk_ρR=ρRnchunk,
        JuliaTask_nchunk_int=Threads.nthreads(),
        BLAS_nthread_epsilon=Threads.nthreads(),
    )

    np = count(ρk.period)
    if isnothing(coulomb)
        if np == 3
            coulomb = PlaneWaveCoulomb(:bulk, ρk)
        elseif np == 2
            coulomb = PlaneWaveCoulomb(:slab, ρk)
        else
            error("You need to specify `PlaneWaveCoulomb` for 1D or 0D!")
        end
    end
    if isnothing(multipolecorrection)
        if np == 3
            println("No multipole correction for 3D, please make sure.")
        end
    else
        if np ≠ 3
            multipolecorrection = nothing
            println("No multipole correction for low-dimensional case, don't use it.")
        end
    end

    BLAS_nt = BLAS.get_num_threads()
    FFTW_nt = FFTW.get_num_threads()
    # 一般使用单线程FFT；单个FFT任务过大可以尝试手动进行1DFFT，以在Julia层级进行更精细的控制。
    FFTW.set_num_threads(1)

    println("[calculate interaction]")
    println("Running on 1 processors.")
    println("nthreads = ", Threads.nthreads())
    _report_rhok(ρk)
    println("plane-wave coulomb: ", string(nameof(typeof(coulomb))))
    println("multipole correction: ", string(nameof(typeof(multipolecorrection))))
    flush(stdout)

    int(Val(np), ρk, iqidx; coulomb, multipolecorrection, outfolder, fileformat, interaction_para)

    println("[calculate interaction end]")

    redirect_stdout(original_stdout)
    close(logio)

    BLAS.set_num_threads(BLAS_nt)
    FFTW.set_num_threads(FFTW_nt)

    return nothing
end
function int(epsilon::DielectricPWBasis, ρk::AbstractRhok;
    ρRnchunk::Integer=Threads.nthreads(), ρRbuffer::Symbol=:auto,
    iqidx=:default, fileformat::Symbol=:jld2,
    outfolder::AbstractString="./", logfile::AbstractString=joinpath(outfolder, "./WannierInterAction.log"),
    interaction::AbstractVector{InterAction}=[DIRECT, EXCHANGE, PAIR],
    coulomb::Union{PlaneWaveCoulomb,Nothing}=nothing,
    multipolecorrection::Union{AbstractMultipoleCorrection,Nothing}=nothing,
    wingaverage::Bool=true,
    grid_atol::Real=1e-6, exchange_range=(3, 3, 3),
    pair=:default, pair_uplo=:all, pair_range=:default, pair_atol::Real=sqrt(eps(Float64)),
)

    mkpath(dirname(logfile))
    logio = get_logio(logfile)
    original_stdout = stdout
    redirect_stdout(logio)

    # used for debugging
    if iqidx == :default
        iqidx = collect(1:ρk.nq)
    else
        iqidx = Int.(collect(iqidx))
        if any(i -> i < 1 || i > ρk.nq, iqidx)
            throw(ArgumentError("iqidx must be between 1 and $(ρk.nq)"))
        end
    end

    interaction = Set(interaction)
    @assert length(interaction) > 0 "Empty `interaction`."

    # check some input parameters that cannot be automatically pre-checked by the program.
    @assert fileformat ∈ [:jld2, :txt] "Wrong `fileformat`!"

    interaction_para = (
        interaction=interaction,
        wingaverage=wingaverage,
        grid_atol=grid_atol,
        exchange_range=exchange_range,
        pair=pair,
        pair_uplo=pair_uplo,
        pair_range=pair_range,
        pair_atol=pair_atol,
        ρR_buffer_level=ρRbuffer,
        JuliaTask_nchunk_ρR=ρRnchunk,
        JuliaTask_nchunk_int=Threads.nthreads(),
        BLAS_nthread_epsilon=Threads.nthreads(),
    )

    np = count(ρk.period)
    if isnothing(coulomb)
        if np == 3
            coulomb = PlaneWaveCoulomb(:bulk, ρk)
        elseif np == 2
            coulomb = PlaneWaveCoulomb(:slab, ρk)
        else
            error("You need to specify `PlaneWaveCoulomb` for 1D or 0D!")
        end
    end
    if isnothing(multipolecorrection)
        if np == 3
            println("No multipole correction for 3D, please make sure.")
        end
    else
        if np ≠ 3
            multipolecorrection = nothing
            println("No multipole correction for low-dimensional case, don't use it.")
        end
    end

    BLAS_nt = BLAS.get_num_threads()
    FFTW_nt = FFTW.get_num_threads()
    # 一般使用单线程FFT；单个FFT任务过大可以尝试手动进行1DFFT，以在Julia层级进行更精细的控制。
    FFTW.set_num_threads(1)

    println("[calculate interaction]")
    println("Running on 1 processors.")
    println("nthreads = ", Threads.nthreads())
    _report_rhok(ρk)
    println("plane-wave coulomb: ", string(nameof(typeof(coulomb))))
    println("multipole correction: ", string(nameof(typeof(multipolecorrection))))
    flush(stdout)

    int(Val(np), epsilon, ρk, iqidx; coulomb, multipolecorrection, outfolder, fileformat, interaction_para)

    println("[calculate interaction end]")

    redirect_stdout(original_stdout)
    close(logio)

    BLAS.set_num_threads(BLAS_nt)
    FFTW.set_num_threads(FFTW_nt)

    return nothing
end
