import shore


server = shore.RemoteServerManager(load='../test.server.pkl')

if __name__ == '__main__':
    # part about the light
    photon = shore.Photon(interaction='dipole', polarization=[0, 0, 1], energy=dict(element='F', edge='K'))
    beam = shore.Light(photons=[photon])

    # part about the material
    lif = shore.AtomicStructure('data/structures/lif.cif')

    material = shore.Matter(
        structure=lif,
        ground_state=dict(bands=10, kpoints=-1, ecut=-1),
        screening=dict(bands=10, kpoints=-1),
        bse=dict(bands=10, kpoints=-1),
        xas=dict(broad=-1, range="1000 -10 40"),
        prefix='srun -n 16'
    )
    

    #creating input
    ocean = shore.Input(name='lio', matter=material, light=beam)
    

    #creating the calculation instance 
    xas = shore.Calculation(input=ocean, server=server)


    #running the calculation
    xas.run()

