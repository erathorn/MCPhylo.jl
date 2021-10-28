using ReverseDiff: GradientTape, GradientConfig, gradient, gradient!, compile, DiffResults

#########
# setup #
#########

# some objective function to work with
f(a, b) = sum(a' * b + a * b')

# pre-record a GradientTape for `f` using inputs of shape 100x100 with Float64 elements
const f_tape = GradientTape(f, (rand(100, 100), rand(100, 100)))

# compile `f_tape` into a more optimized representation
const compiled_f_tape = compile(f_tape)

# some inputs and work buffers to play around with
a, b = rand(100, 100), rand(100, 100)
inputs = (a, b)
results = (similar(a), similar(b))
all_results = map(DiffResults.GradientResult, results)
cfg = GradientConfig(inputs)

####################
# taking gradients #
####################

# with pre-recorded/compiled tapes (generated in the setup above) #
#-----------------------------------------------------------------#

# this should be the fastest method, and non-allocating
gradient!(results, compiled_f_tape, inputs)

# the same as the above, but in addition to calculating the gradients, the value `f(a, b)`
# is loaded into the the provided `DiffResult` instances (see DiffResults.jl documentation).
gradient!(all_results, compiled_f_tape, inputs)

# this should be the second fastest method, and also non-allocating
gradient!(results, f_tape, inputs)

# you can also make your own function if you want to abstract away the tape
∇f!(results, inputs) = gradient!(results, compiled_f_tape, inputs)

# with a pre-allocated GradientConfig #
#-------------------------------------#
# these methods are more flexible than a pre-recorded tape, but can be
# wasteful since the tape will be re-recorded for every call.

gradient!(results, f, inputs, cfg)

gradient(f, inputs, cfg)

# without a pre-allocated GradientConfig #
#----------------------------------------#
# convenient, but pretty wasteful since it has to allocate the GradientConfig itself

gradient!(results, f, inputs)

gradient(f, inputs)



[ 2.0450003722816574 -0.4036882171184115 -0.29246840503742927 0.5394489737681402 0.2635587384209584;
 -0.4036882171184115 2.1945387508796306 0.33001618120050635 -0.6769908583435065 -0.35865598799698073;
 -0.29246840503742927 0.33001618120050635 1.6405268101691983 -0.561004898578046 -0.28755129934965756;
  0.5394489737681402 -0.6769908583435065 -0.561004898578046 4.400834882333025 0.8934881192114935;
  0.2635587384209584 -0.35865598799698073 -0.2875512993496576 0.8934881192114936 1.8709900901568826]

  [0.14186354708749327 0.04803891080994312 0.01857073967875705 0.08967903942348894 -0.13727180581043386;
   0.04803891080994312 0.15069653217754292 0.08668300648124458 -0.1111310230742415 -0.12469643909821908;
   0.01857073967875705 0.08668300648124458 0.3028455106352376 -0.17410349472045375 -0.14769542287057527;
   0.08967903942348894 -0.1111310230742415 -0.17410349472045375 0.6004764415185736 -0.025912730971656226;
   -0.13727180581043386 -0.12469643909821908 -0.14769542287057527 -0.025912730971656226 0.32612112320314296]