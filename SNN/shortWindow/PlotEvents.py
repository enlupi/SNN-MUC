import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({'figure.max_open_warning': 0})
from Params import *
from EventData import *


pos_z = [ np.nan ] + [ZCELL*(-(NLAYERS-1-(NLAYERS%2))/2 + i) for i in range(NLAYERS)]
xmin = -(1+(NWIRES-1)//2)*XCELL
pos_x = [ np.nan ] + [xmin + (1+s)*0.5*XCELL for s in is_shifted_right]


def plot_muon(muon_hits, bx0, x0=None, m=None, t0_scint=None, save_me=False):
    """
    Plot a list of muon hits
    Args:
        - muon_hits: list of dicts containing at least the following
                     information:
                        - hit layer
                        - hit wire (in the macrocell)
                        - hit bx
                        - hit tdc
                        - hit side (-1,0,1)
        - bx0: time pedestal for this macrocell
        - x0: local intercept of the track
        - m: tan_psi of the track
    """

    fig, ax = plt.subplots(1, 1, figsize=(XCELL*(NWIRES+1)/10,NLAYERS*ZCELL/10))

    def the_grid(the_plot):
        for ilay in range(1, NLAYERS+1):
            for iwire in range(1, NWIRES+1):
                the_plot.plot(pos_x[ilay]+XCELL*(iwire-1),pos_z[ilay],
                              marker='.',
                              markersize=2,
                              color='grey',
                              zorder=-10)
                the_plot.add_patch(
                    Rectangle((pos_x[ilay]-0.5*XCELL+XCELL*(iwire-1),
                                   pos_z[ilay]-ZCELL*0.5),
                                  XCELL,
                                  ZCELL,
                                 edgecolor='lightgrey',
                                 facecolor='None',
                             zorder=-10)
                )

    def the_planes(the_plot):
        for ilay in  range(1, NLAYERS+1):
            the_plot.add_patch(
                Rectangle((pos_x[ilay]-0.5*XCELL,
                           pos_z[ilay]-ZCELL*0.5-PLANE_WIDTH/2.),
                          NWIRES*XCELL,
                          PLANE_WIDTH,
                          facecolor='darkgrey',
                          edgecolor='None',
                          zorder=-15)
            )
            the_plot.add_patch(
                Rectangle((pos_x[ilay]-0.5*XCELL,
                           ZCELL+pos_z[ilay]-ZCELL*0.5-PLANE_WIDTH/2.),
                          NWIRES*XCELL,
                          PLANE_WIDTH,
                          facecolor='darkgrey',
                          edgecolor='None',
                          zorder=-15)
            )


    def the_ibeams(the_plot):
        for ilay in  range(1, NLAYERS+1):
            for iwire in range(1, NWIRES+1):
                the_plot.add_patch(
                    Rectangle((pos_x[ilay]-0.5*XCELL+XCELL*(iwire-1)-IBEAM_WIDTH/2,
                               pos_z[ilay]-ZCELL*0.5-PLANE_WIDTH/2.),
                              IBEAM_WIDTH,
                              ZCELL+PLANE_WIDTH,
                              facecolor='darkgrey',
                              edgecolor='None',
                              zorder=-13)
                )
                the_plot.add_patch(
                    Rectangle((XCELL+pos_x[ilay]-0.5*XCELL+XCELL*(iwire-1)-IBEAM_WIDTH/2,
                               pos_z[ilay]-ZCELL*0.5-PLANE_WIDTH/2.),
                              IBEAM_WIDTH,
                              ZCELL+PLANE_WIDTH,
                              facecolor='darkgrey',
                              edgecolor='None',
                              zorder=-13)
                )

                the_plot.add_patch(
                    Rectangle((XCELL+pos_x[ilay]-0.5*XCELL+XCELL*(iwire-1)-IBEAM_WING/2,
                               pos_z[ilay]-ZCELL*0.5+PLANE_WIDTH/2.),
                              IBEAM_WING,
                              IBEAM_WIDTH,
                              facecolor='darkgrey',
                              edgecolor='None',
                              zorder=-13)
                )
                the_plot.add_patch(
                    Rectangle((XCELL+pos_x[ilay]-0.5*XCELL+XCELL*(iwire-1)-IBEAM_WING/2,
                                  pos_z[ilay]+ZCELL*0.5-PLANE_WIDTH/2.-IBEAM_WIDTH),
                                  IBEAM_WING,
                                  IBEAM_WIDTH,
                                  facecolor='darkgrey',
                                  edgecolor='None',
                             zorder=-13)
                )

                the_plot.add_patch(
                    Rectangle((XCELL+pos_x[ilay]-0.5*XCELL+XCELL*(iwire-2)-IBEAM_WING/2,
                                  pos_z[ilay]-ZCELL*0.5+PLANE_WIDTH/2.),
                                  IBEAM_WING,
                                  IBEAM_WIDTH,
                                  facecolor='darkgrey',
                                  edgecolor='None',
                             zorder=-13)
                )
                the_plot.add_patch(
                    Rectangle((XCELL+pos_x[ilay]-0.5*XCELL+XCELL*(iwire-2)-IBEAM_WING/2,
                                  pos_z[ilay]+ZCELL*0.5-PLANE_WIDTH/2.-IBEAM_WIDTH),
                                  IBEAM_WING,
                                  IBEAM_WIDTH,
                                  facecolor='darkgrey',
                                  edgecolor='None',
                             zorder=-13)
                )



    the_grid(ax)
    the_planes(ax)
    the_ibeams(ax)

    for hit in muon_hits:
        layer, wire, bx, tdc =  hit['layer'], hit['wire_num'], hit['bx'], hit['tdc']
        plt.scatter(pos_x[layer] + (wire-1)*XCELL, pos_z[layer], color='red')

        # x left and right
        wire_pos = pos_x[layer] + (wire-1)*XCELL
        wire_z   = pos_z[layer]

        tdrift = (bx + tdc/30 - bx0)*DURATION['bx']
        x_drift = tdrift * VDRIFT
        """print("z: {}".format(wire_z))
        print("x-: {}".format(wire_pos - x_drift))
        print("x+: {}".format(wire_pos + x_drift))"""
        plt.scatter(wire_pos - x_drift, wire_z, color='green', marker="x")
        plt.scatter(wire_pos + x_drift, wire_z, color='green', marker="x")

        if t0_scint != None:
            tdrift = (bx + tdc/30 - t0_scint)*DURATION['bx']
            x_drift = tdrift * VDRIFT
            plt.scatter(wire_pos - x_drift, wire_z,  s=80, facecolors='none', edgecolors='b')
            plt.scatter(wire_pos + x_drift, wire_z,  s=80, facecolors='none', edgecolors='b')


    if (x0!=None) and (m!=None):
        bottom_x = x0 + m*(pos_z[1]-0.5*ZCELL)
        top_x    = x0 + m*(pos_z[NLAYERS]+0.5*ZCELL)

        plt.plot((bottom_x,top_x),(pos_z[1]-0.5*ZCELL,pos_z[NLAYERS]+0.5*ZCELL),'b-')


def plot_event(event):
    muon_hits = numpy_to_hits(event)
    m, x0 = np.tan(muon_hits[0]['psi']), muon_hits[0]['x0']
    plot_muon(muon_hits=muon_hits, bx0=500, m=m, x0=x0)