using StructTypes
struct Person
    id::Int
    name::String
end

StructTypes.StructType(::Type{Person}) = StructTypes.Struct()
struct PersonWrapper
    person::Person
    second_person::Int64
end
StructTypes.StructType(::Type{PersonWrapper}) = StructTypes.CustomStruct()
StructTypes.lower(x::PersonWrapper) = [x.person, x.second_person]
StructTypes.lowertype(::Type{PersonWrapper}) = Vector{Any}
StructTypes.construct(::Type{PersonWrapper}, x::Vector{Any}) = PersonWrapper(x...)

using JSON3
pw = PersonWrapper(Person(1, "blub"), 2)
s = JSON3.write(pw)
t = JSON3.read(s, PersonWrapper)
