### To test additions to the package-

1) In the "doc" folder, either create a subfolder, OR create a file in the relevant pre-existing subfolder.
	(to provide tests for the implementation of a new parser, add your tests to the "parsers" subfolder, for example)

2) Place unit test file and all relevant helper files in your chosen subfolder.

3) Navigate to the "test" folder.

4) In runtests.jl, format your tests as follows:
First, add the filename of your unit test file to the const that corresponds to the subfolder it is in;

for example:
```Julia
		const parsertests = [
			"newick"
			]
```
where "newick" is the name of a file containing tests in the "parsers" subfolder. If a new subfolder has been created, it will be necessary to add a corresponding const.

If your test file was added to a preexisting subfolder, no further steps are necessary. If a new subfolder was created, add an additional @testset clause with the following format:

```Julia
    @testset "NAMEOFTESTS" begin
        for t in NAMEOFCONST
        @everywhere Random.seed!(123)
		@runtest "../doc/NAMEOFSUBFOLDER/" t
	  end
	  end
