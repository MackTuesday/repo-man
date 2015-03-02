/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yinyangfire;

import java.awt.Canvas;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.Graphics2D;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsEnvironment;
import java.awt.image.BufferStrategy;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.Toolkit;
import java.awt.Transparency;

import java.util.Random;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

public class YinYangFire extends Thread {
    private boolean isRunning = true;
    private Canvas canvas;
    private BufferStrategy strategy;
    private BufferedImage background;
    private Graphics2D backgroundGraphics;
    private Graphics2D graphics;
    private JFrame frame;
    private final int width = 512;
    private final int height = 512;
    private GraphicsConfiguration config =
    		GraphicsEnvironment.getLocalGraphicsEnvironment()
                                   .getDefaultScreenDevice()
                                   .getDefaultConfiguration();
    private int frameCount = 0;
    private int[][][] cellStates;
    private int currentCellArray = 0;
    private final int numStates = 64;

    // Setup
    public YinYangFire() {
    	// JFrame
    	frame = new JFrame();
    	frame.addWindowListener(new FrameClose());
    	frame.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
    	frame.setSize(width, height);
    	frame.setVisible(true);

    	// Canvas
    	canvas = new Canvas(config);
    	canvas.setSize(width, height);
    	frame.add(canvas, 0);

    	// Background & Buffer
    	background = create(width, height, false);
    	canvas.createBufferStrategy(2);
    	do {
            strategy = canvas.getBufferStrategy();
    	} while (strategy == null);
        
        cellStates = new int[width][height][2];
        Random randGen = new Random();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
               cellStates[x][y][currentCellArray] = randGen.nextInt(numStates);
            }
        }
        
    	start();
    }

    // create a hardware accelerated image
    public final BufferedImage create(final int width, final int height,
                                      final boolean alpha) {
    	return config.createCompatibleImage(width, height,
                alpha ? Transparency.TRANSLUCENT : Transparency.OPAQUE);
    }

    public ColorModel getColorModel ()
	      {
        return ColorModel.getRGBdefault();
    }

    private class FrameClose extends WindowAdapter {
    	@Override
    	public void windowClosing(final WindowEvent e) {
            isRunning = false;
    	}
    }

    // Screen and buffer stuff
    private Graphics2D getBuffer() {
    	if (graphics == null) {
            try {
                graphics = (Graphics2D) strategy.getDrawGraphics();
            } catch (IllegalStateException e) {
                return null;
            }
    	}
    	return graphics;
    }

    private boolean updateScreen() {
    	graphics.dispose();
    	graphics = null;
    	try {
            strategy.show();
            Toolkit.getDefaultToolkit().sync();
            return (!strategy.contentsLost());
        } catch (NullPointerException e) {
            return true;
        } catch (IllegalStateException e) {
            return true;
    	}
    }

    @Override
    public void run() {
    	backgroundGraphics = (Graphics2D) background.getGraphics();
    	long fpsWait = (long) (1.0 / 30 * 1000);
    	main: while (isRunning) {
            long renderStart = System.nanoTime();
            updateGame();

            // Update Graphics
            do {
                Graphics2D bg = getBuffer();
                if (!isRunning) {
                    break main;
                }
                renderGame(backgroundGraphics); // this is your draw method
                bg.drawImage(background, 0, 0, width, height,
                                         0, 0, width, height, null);
                bg.dispose();
            } while (!updateScreen());

            // Better do some FPS limiting here
            long renderTime = (System.nanoTime() - renderStart) / 1000000;
            /*try {
                Thread.sleep(Math.max(0, fpsWait - renderTime));
            } catch (InterruptedException e) {
                Thread.interrupted();
                break;
            }*/
            
            frameCount++;
    	}
        
    	frame.dispose();
    }

    public void updateGame() {
    	int nextCellArray = 1 - currentCellArray;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int lt = (x-1+width)  % width;
                int rt = (x+1+width)  % width;
                int up = (y-1+height) % height;
                int dn = (y+1+height) % height;
                int stateSum = cellStates[lt][up][currentCellArray] +
                               cellStates[ x][up][currentCellArray] +
                               cellStates[rt][up][currentCellArray] +
                               cellStates[lt][ y][currentCellArray] +
                               cellStates[rt][ y][currentCellArray] +
                               cellStates[lt][dn][currentCellArray] +
                               cellStates[ x][dn][currentCellArray] +
                               cellStates[rt][dn][currentCellArray];
                int me = cellStates[ x][ y][currentCellArray];
                if (me * 8 + 2 >= stateSum) {
                    me--;
                    if (me < 0) {
                        me += numStates;
                    }
                } else {
                    me++;
                }
                
                cellStates[x][y][nextCellArray] = me;
            }
        }
        
        currentCellArray = nextCellArray;
    }

    public void renderGame(Graphics2D g) {
        int rgbArray[] = new int[width*height];
//        for (int i = 0; i < width*height; i++) {
//            rgbArray[i] = (i + frameCount*width*16);
//        }
        int arrayIndex = 0;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                rgbArray[arrayIndex] = cellStates[x][y][currentCellArray] * 4;
                arrayIndex++;
            }
        }
        background.setRGB(0, 0, width, height, rgbArray, 0, width);
    }

    public static void main(final String args[]) {
    	new YinYangFire();
    }
}
