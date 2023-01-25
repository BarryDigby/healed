a = [ [ id:foo, tab:a ], file(main.nf)]
     

b = [ [ id:foo, tab:b ], file(main.nf)]


c = channel.of(a)
d = channel.of(b)

z = c.merge(d).view()
