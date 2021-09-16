struct FileSyntaxError <: Exception
    info::String
  end
  
  Base.showerror(io::IO, e::FileSyntaxError) = print(io, "FileSyntaxError: ", e.info)  