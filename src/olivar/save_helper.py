import logging
logger = logging.getLogger('main')

import json
import os
import pickle
import pandas as pd
import numpy as np

from Bio import SeqIO

import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots

DESIGNEXT = '.olvd' # extension for Olivar design file


def save(design_out, out_path: str):
    '''
    Load from a previous Olivar design object (in tiling()) or file (.olvd) and save output files.
    Input:
        design_out: Output of tiling(), or load the olvd file with pickle.
        out_path: Output directory.
    '''
    # load data
    if type(design_out) is str:
        if os.path.isfile(design_out) and design_out.endswith(DESIGNEXT):
            logger.info(f'Loading Olivar design from {design_out}...')
            with open(design_out,  'rb') as f:
                design_out = pickle.load(f)
                logger.info(f'Successfully loaded Olivar design.')
        else:
            raise FileNotFoundError(f'Olivar design (.olvd) file "{design_out}" not found or is invalid.')
    
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        
    config = design_out['config']
    lc = design_out['learning_curve']
    df = design_out['df']
    all_ref_info = design_out['all_ref_info']

    # column name mapper for backward compatibility
    df.rename(columns={
        'amp_id': 'amplicon_id', 
        'amp': 'amplicon'
    }, inplace=True)

    # save configurations
    save_path = os.path.join(out_path, '%s.json' % config['title'])
    with open(save_path, 'w') as f:
        json.dump(config, f, indent=4)
    logger.info(f'Configurations saved as {save_path}')

    # save sequences and coordinates as csv
    save_path = os.path.join(out_path, '%s.csv' % config['title'])
    df.to_csv(save_path, index=False)
    logger.info(f'Primer pools saved as {save_path}')

    # save ARTIC output format
    try:
        art_df = design_out['art_df']
        save_path = os.path.join(out_path, '%s.primer.bed' % config['title'])
        art_df.sort_values(by='start', axis=0, ascending=True, inplace=True)
        art_df.to_csv(save_path, sep='\t', header=False, index=False)
        logger.info(f'Primer pools (ARTIC format) saved as {save_path}')
    except KeyError:
        logger.info('ARTIC/PrimalScheme format is not included in older versions of Olivar, skipped.')

    #------------- SADDLE Loss -------------#
    fig = make_subplots(
        rows=1, cols=2, 
        subplot_titles=('Primer dimer optimization (pool-1)', 
            'Primer dimer optimization (pool-2)')
    )

    fig.add_trace(
        go.Scatter(
            y=lc[0], 
            line=dict(color='#1f77b4'), 
            showlegend=False
        ), 
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(
            y=lc[1], 
            line=dict(color='#1f77b4'), 
            showlegend=False
        ), 
        row=1, col=2
    )

    fig.update_xaxes(title_text='iterations', row=1, col=1)
    fig.update_xaxes(title_text='iterations', row=1, col=2)

    fig.update_yaxes(title_text='SADDLE Loss', row=1, col=1)
    fig.update_yaxes(title_text='SADDLE Loss', row=1, col=2)

    # save html figure
    save_path = os.path.join(out_path, f"{config['title']}_SADDLE_Loss.html")
    with open(save_path, 'w') as f:
        f.write(plotly.io.to_html(fig))
    logger.info(f'SADDLE optimization plot saved as {save_path}')
    #------------- SADDLE Loss -------------#

    for ref_name, ref_info in all_ref_info.items():
        logger.info(f'Saving output files and figures for {ref_name}...')

        risk_arr = ref_info['risk_arr']
        gc_arr = ref_info['gc_arr']
        comp_arr = ref_info['comp_arr']
        hits_arr = ref_info['hits_arr']
        var_arr = ref_info['var_arr']
        seq_record = ref_info['seq_record']
        all_loss = ref_info['all_loss']

        # save reference sequence
        save_path = os.path.join(out_path, '%s_ref.fasta' % ref_name)
        with open(save_path, 'w') as f:
            SeqIO.write([seq_record], f, 'fasta')
        logger.info(f'Reference sequence saved as {save_path}')

        # save risk array
        risk = pd.DataFrame({
            'position': range(1, len(risk_arr)+1), 
            'base': list(seq_record.seq), 
            'extreme GC': gc_arr, 
            'low complexity': comp_arr, 
            'non-specificity': hits_arr, 
            'variations': var_arr, 
            'risk': risk_arr
        })
        save_path = os.path.join(out_path, '%s_risk.csv' % ref_name)
        risk.to_csv(save_path, index=False)
        logger.info(f'Risk scores saved as {save_path}')

        #------------- PDR Loss -------------#
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                y=all_loss, 
                line=dict(color='#1f77b4'), 
                showlegend=False
            )
        )
        fig.update_layout(
            title=f'Optimization of primer design regions (PDRs) for {ref_name}', 
            xaxis_title='iterations (sorted)', 
        )
        fig.update_yaxes(title_text='Loss', type='log')
        # save html figure
        save_path = os.path.join(out_path, f'{ref_name}_PDR_Loss.html')
        with open(save_path, 'w') as f:
            f.write(plotly.io.to_html(fig))
        logger.info(f'PDR optimization plot saved as {save_path}')
        #------------- PDR Loss -------------#

        #------------- risk array and primers -------------#
        # create figure
        fig = go.Figure()

        # set figure scale
        base_offset = (0.5*config['w_egc'] + 0.5*config['w_lc'] + config['w_ns'])/6
        base_offset *= 3

        # plot risk array
        r = np.arange(len(risk_arr))
        fig.add_trace(
            go.Scatter(
                x=r, y=gc_arr,
                hoverinfo='skip',
                mode='lines',
                line=dict(width=0, color='#1f77b4'),
                name='extreme GC',
                stackgroup='one' # define stack group
            )
        )
        fig.add_trace(
            go.Scatter(
                x=r, y=comp_arr,
                hoverinfo='skip',
                mode='lines',
                line=dict(width=0, color='#ff7f0e'),
                name='low complexity',
                stackgroup='one'
            )
        )
        fig.add_trace(
            go.Scatter(
                x=r, y=hits_arr,
                hoverinfo='skip',
                mode='lines',
                line=dict(width=0, color='#2ca02c'),
                name='non-specificity',
                stackgroup='one'
            )
        )
        fig.add_trace(
            go.Scatter(
                x=r, y=var_arr,
                hoverinfo='skip',
                mode='lines',
                line=dict(width=0, color='#d62728'),
                name='variations',
                stackgroup='one'
            )
        )

        # plot primers
        fp_rp_diff = 0.1 # distance between fP and rP
        head_dx = 7 # arrow head length
        head_dy = base_offset*7/90 # arrow head height
        primer_x = [] # x coords for primer plot
        primer_y = [] # y coords for primer plot
        hover_text = []

        # pool 1
        fp_offset = 2 * base_offset
        rp_offset = (2-fp_rp_diff) * base_offset
        for i, row in df[(df['reference']==ref_name) & (df['pool']==1)].iterrows():
            fp_start = row['start']-1
            fp_stop = row['insert_start']-1
            rp_start = row['insert_end']
            rp_stop = row['end']
            # [fp_start, fp_stop] and [fp_offset, fp_offset] is the body of fP
            # [fp_stop-head_dx, fp_stop] and [fp_offset+head_dy, fp_offset] is the head of fP
            primer_x.extend([fp_start, fp_stop, None, fp_stop-head_dx, fp_stop, None, 
                rp_start, rp_stop, None, rp_start, rp_start+head_dx, None])
            primer_y.extend([fp_offset, fp_offset, None, fp_offset+head_dy, fp_offset, None, 
                rp_offset, rp_offset, None, rp_offset, rp_offset-head_dy, None])
            hover_text.extend(['pool-1 %s fP: %d - %d' % (row['amplicon_id'], row['start'], fp_stop)]*6 + \
                ['pool-1 %s rP: %d - %d' % (row['amplicon_id'], rp_start+1, rp_stop)]*6)
        
        # pool 2
        fp_offset = base_offset
        rp_offset = (1-fp_rp_diff) * base_offset
        for i, row in df[(df['reference']==ref_name) & (df['pool']==2)].iterrows():
            fp_start = row['start']-1
            fp_stop = row['insert_start']-1
            rp_start = row['insert_end']
            rp_stop = row['end']
            # [fp_start, fp_stop] and [fp_offset, fp_offset] is the body of fP
            # [fp_stop-head_dx, fp_stop] and [fp_offset+head_dy, fp_offset] is the head of fP
            primer_x.extend([fp_start, fp_stop, None, fp_stop-head_dx, fp_stop, None, 
                rp_start, rp_stop, None, rp_start, rp_start+head_dx, None])
            primer_y.extend([fp_offset, fp_offset, None, fp_offset+head_dy, fp_offset, None, 
                rp_offset, rp_offset, None, rp_offset, rp_offset-head_dy, None])
            hover_text.extend(['pool-2 %s fP: %d - %d' % (row['amplicon_id'], row['start'], fp_stop)]*6 + \
                ['pool-2 %s rP: %d - %d' % (row['amplicon_id'], rp_start+1, rp_stop)]*6)
        
        # plot all primers
        fig.add_trace(
            go.Scatter(
                x=primer_x, y=primer_y, 
                mode='lines', # connect points with lines
                connectgaps=False, # gaps (np.nan or None) in the provided data arrays are not connected
                line=dict(color='rgb(0,0,0)'), # black
                showlegend=False, 
                hovertemplate='%{text}<extra></extra>', # <extra></extra> hides the "trace 4"
                text=hover_text
            )
        )

        # title and axis
        fig.update_layout(
            hoverlabel=dict(
                bgcolor='white'
            ), 
            # title
            title={
                'text': f"{ref_name} (risk components are stacked together)",
                'x': 0.98, 
                'font': {
                    'size': 18,
                }
            }, 

            # axis
            xaxis=dict(
                tickformat='%d', 
                tickfont=dict(
                    size=14, 
                ), 
                showgrid=False, 
                rangeslider=dict(
                    visible=True
                )
            ), 
            yaxis=dict(
                range=[0, base_offset*6], 
                tickfont=dict(
                    size=14, 
                ), 
                showgrid=False
            ), 

            # axis title
            xaxis_title='position', 
            yaxis_title='risk', 

            # global font
            font=dict(
                size=16,
            )
        )

        # save html figure
        save_path = os.path.join(out_path, '%s.html' % ref_name)
        with open(save_path, 'w') as f:
            f.write(plotly.io.to_html(fig))
        logger.info(f'Risk and primer viewer saved as {save_path}')
        #------------- risk array and primers -------------#
