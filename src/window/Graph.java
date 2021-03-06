package window;

import solution.Data;
import solution.Solution;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * @author misha
 */
public class Graph extends javax.swing.JFrame {
    private int countOfGraph = 0;
    private final XYSeriesCollection xysc = new XYSeriesCollection( );
    private boolean isFirstGraph = true;

    public void addGraph( 
            double[] X, double[] Y, String name) {
        XYSeries series = new XYSeries( name);
        for ( int i = 0; i < X.length; i++)
            series.add( X[i], Y[i]);
        xysc.addSeries( series);
    }

    public void clear( ) {
        xysc.removeAllSeries( );
        countOfGraph = 0;
        drawPanel.validate( );
        drawPanel.updateUI( );
    }

    public void building( String nameMethods) {
        JFreeChart chart = ChartFactory.createXYLineChart( 
                nameMethods, "backslashu03B8", "u",
                xysc,
                PlotOrientation.VERTICAL,
                true, true, true);
        drawPanel.setLayout( new java.awt.BorderLayout( ));
        drawPanel.add( new ChartPanel( chart));
        drawPanel.validate( );
        drawPanel.updateUI( );
    }

    /**
     * Creates new form Graph
     */
    public Graph( ) {
        initComponents( );
        startButton.addActionListener( new ActionListener( ) {
            @Override
            public void actionPerformed( ActionEvent actionEvent) {
                startButtonActionPerformed( actionEvent);
            }
        });
        clearButton.addActionListener( new ActionListener( ) {
            @Override
            public void actionPerformed( ActionEvent actionEvent) {
                clearButtonActionPerformed( actionEvent);
            }
        });
    }
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings( "unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">
    private void initComponents( ) {
        cTextField = new javax.swing.JTextField( );
        cFormattedTextField = new javax.swing.JFormattedTextField( );
        kFormattedTextField = new javax.swing.JFormattedTextField( );
        rFormattedTextField = new javax.swing.JFormattedTextField( );
        betaFormattedTextField = new javax.swing.JFormattedTextField( );
        epsFormattedTextField = new javax.swing.JFormattedTextField( );
        startButton = new javax.swing.JButton( );
        rTextField = new javax.swing.JTextField( );
        KTextField = new javax.swing.JTextField( );
        epsTextField = new javax.swing.JTextField( );
        betaTextField = new javax.swing.JTextField( );
        drawPanel = new javax.swing.JPanel( );
        tTextField = new javax.swing.JTextField( );
        tFormattedTextField = new javax.swing.JFormattedTextField( );
        clearButton = new javax.swing.JButton( );

        setDefaultCloseOperation( javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setName( "basicFrame"); // NOI18N

        cTextField.setEditable( false);
        cTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        cTextField.setText( "c");

        cFormattedTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        cFormattedTextField.addActionListener( new java.awt.event.ActionListener( ) {
            public void actionPerformed( java.awt.event.ActionEvent evt) {
                cFormattedTextFieldActionPerformed( evt);
            }
        });
        cFormattedTextField.setText( "2.0");

        kFormattedTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        kFormattedTextField.addActionListener( new java.awt.event.ActionListener( ) {
            public void actionPerformed( java.awt.event.ActionEvent evt) {
                kFormattedTextFieldActionPerformed( evt);
            }
        });
        kFormattedTextField.setText( "0.01");

        rFormattedTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        rFormattedTextField.addActionListener( new java.awt.event.ActionListener( ) {
            public void actionPerformed( java.awt.event.ActionEvent evt) {
                rFormattedTextFieldActionPerformed( evt);
            }
        });
        rFormattedTextField.setText( "2.0");

        betaFormattedTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        betaFormattedTextField.addActionListener( new java.awt.event.ActionListener( ) {
            public void actionPerformed( java.awt.event.ActionEvent evt) {
                betaFormattedTextFieldActionPerformed( evt);
            }
        });
        betaFormattedTextField.setText( "0.1");

        epsFormattedTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        epsFormattedTextField.addActionListener( new java.awt.event.ActionListener( ) {
            public void actionPerformed( java.awt.event.ActionEvent evt) {
                epsFormattedTextFieldActionPerformed( evt);
            }
        });
        epsFormattedTextField.setText( "1e-7");

        startButton.setText( "start");

        clearButton.setText( "clear");

        rTextField.setEditable( false);
        rTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        rTextField.setText( "R");

        KTextField.setEditable( false);
        KTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        KTextField.setText( "K");

        epsTextField.setEditable( false);
        epsTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        epsTextField.setText( "backslashu03B5");

        betaTextField.setEditable( false);
        betaTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        betaTextField.setText( "backslashu03B2");

        tTextField.setEditable( false);
        tTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        tTextField.setText( "t");

        tFormattedTextField.setHorizontalAlignment( javax.swing.JTextField.CENTER);
        tFormattedTextField.addActionListener( new java.awt.event.ActionListener( ) {
            public void actionPerformed( java.awt.event.ActionEvent evt) {
                tFormattedTextFieldActionPerformed( evt);
            }
        });
        tFormattedTextField.setText( "0.001");

        drawPanel.setBackground( new java.awt.Color( 255, 255, 255));

        javax.swing.GroupLayout drawPanelLayout = new javax.swing.GroupLayout( drawPanel);
        drawPanel.setLayout( drawPanelLayout);
        drawPanelLayout.setHorizontalGroup( 
              //...
        );
        layout.setVerticalGroup( 
              //...
        );

        pack( );
    }// </editor-fold>

    // ... 
    private void startButtonActionPerformed( java.awt.event.ActionEvent evt) {
        // TODO add your handling code here:
        try {
            Data.beta_doleq( Double.parseDouble( betaFormattedTextField.getText( )));
            Data.c_doleq( Double.parseDouble( cFormattedTextField.getText( )));
            Data.eps_doleq( Double.parseDouble( epsFormattedTextField.getText( )));
            Data.K_doleq( Double.parseDouble( kFormattedTextField.getText( )));
            Data.R_doleq( Double.parseDouble( rFormattedTextField.getText( )));
            Data.t_doleq( Double.parseDouble( tFormattedTextField.getText( )));

            double t = Data.t( );
            scala.Tuple2<double[],double[]> res = Solution.masU( t,100);
            double[] X = res._1( );
            double[] Y = res._2( );

            java.text.DecimalFormat f3 = new java.text.DecimalFormat( "0.000");
            addGraph( X, Y, "u" + ( ( countOfGraph == 0) ? "" : countOfGraph) +
                    "( " + f3.format( t) + "; θ)");
            countOfGraph++;
            if ( isFirstGraph) {
                isFirstGraph = false;
                building( "Graph of function u( t,θ)");
            }
        }  catch ( java.lang.NumberFormatException | java.lang.NullPointerException e){
            javax.swing.JOptionPane.showMessageDialog( this,"Running impossible: incorrect data",
                "Error",javax.swing.JOptionPane.WARNING_MESSAGE);
        }
}
    private void startButtonMouseClicked( java.awt.event.ActionEvent evt) {
        startButtonActionPerformed( evt);
    }

    private void clearButtonActionPerformed( java.awt.event.ActionEvent evt) {
        // TODO add your handling code here:
        clear( );
    }
    /**
     * @param args the command line arguments
     */
    public static void main( String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code ( optional) ">
        /* If Nimbus ( introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html
         */
        try {
            for ( javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels( )) {
                if ( "Nimbus".equals( info.getName( ))) {
                    javax.swing.UIManager.setLookAndFeel( info.getClassName( ));
                    break;
                }
            }
        } catch ( ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger( Graph.class.getName( )).log( java.util.logging.Level.SEVERE, null, ex);
        } catch ( InstantiationException ex) {
            java.util.logging.Logger.getLogger( Graph.class.getName( )).log( java.util.logging.Level.SEVERE, null, ex);
        } catch ( IllegalAccessException ex) {
            java.util.logging.Logger.getLogger( Graph.class.getName( )).log( java.util.logging.Level.SEVERE, null, ex);
        } catch ( javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger( Graph.class.getName( )).log( java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater( new Runnable( ) {
            public void run( ) {
                new Graph( ).setVisible( true);
            }
        });
    }

    // Variables declaration - do not modify
    private javax.swing.JTextField KTextField;
    private javax.swing.JFormattedTextField betaFormattedTextField;
    private javax.swing.JTextField betaTextField;
    private javax.swing.JFormattedTextField cFormattedTextField;
    private javax.swing.JTextField cTextField;
    private javax.swing.JButton clearButton;
    private javax.swing.JPanel drawPanel;
    private javax.swing.JFormattedTextField epsFormattedTextField;
    private javax.swing.JTextField epsTextField;
    private javax.swing.JFormattedTextField kFormattedTextField;
    private javax.swing.JFormattedTextField rFormattedTextField;
    private javax.swing.JTextField rTextField;
    private javax.swing.JButton startButton;
    private javax.swing.JFormattedTextField tFormattedTextField;
    private javax.swing.JTextField tTextField;
    // End of variables declaration
}
