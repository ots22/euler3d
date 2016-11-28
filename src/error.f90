module m_error
contains
  subroutine warn(msg)
    character(len=*), intent(in) :: msg
    write (0,*) "warning: ", msg
  end subroutine warn

  subroutine panic(msg)
    character(len=*), intent(in) :: msg
    write (0,*) msg
    error stop
  end subroutine panic

  subroutine assert(expr,msg)
    logical, intent(in) :: expr
    character(len=*), intent(in) :: msg
    if (.not.expr) call panic("Assert failure: " // msg)
    continue
  end subroutine assert
end module m_error
