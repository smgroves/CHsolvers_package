function test(a; b=[1, 0, 1, 0])
    println(a)
    xleft = b[1]
    xright = b[2]
    println(xleft)
    println(xright)
    xleft, xright, yleft, yright = b
    println(yleft)
    println(yright)
end

test(1, b=[2, 3, 4, 5])
