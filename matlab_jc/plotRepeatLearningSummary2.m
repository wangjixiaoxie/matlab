function plotRepeatLearningSummary2(Exp,probstay,numstay,alls2,indPost)
                  figure;hold on;
                    subplot(6,4,17);hold on;
                  for i=[4 6 10 12]
                      if Exp(i).hitabove==1
                          if ~Exp(i).postlesion
                              for m=1:length(probstay(i,1).rn)
                                  if min([numstay(i,1).rn(m) numstay(i,2).rn(m)])>10
                                      if m>Exp(i).MIN-2 % if it's hit 
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ro','Markersize',5)
                                      else
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ko','Markersize',5) 
                                      end
                                  end
                              end
                          end
                      end
                  end
                  plot([0 1],[0 1],'k')
                  subplot(6,4,18);hold on;
                  plot(alls2([4 6 10 12],1),alls2([4 6 10 12],2),'ro','Markersize',5)
                  plot([0 10],[0 10],'k')

                 subplot(6,4,21);hold on;
                  for i=[4 6 10 12]
                      if Exp(i).hitabove==1
                          if ~Exp(i).postlesion
                              for m=1:length(probstay(i,1).rn)
                                  if min([numstay(i,1).rn(m) numstay(i,2).rn(m)])>10
                                      if m>Exp(i).MIN-2 % if it's hit 
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ro','Markersize',5)
                                      else
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ko','Markersize',5) 
                                      end
                                      if i~=6
                                      plot(probstay(i,1).rn(m),probstay(i,3).rn(m),'bo','Markersize',5) 
                                      end
                                  end
                              end
                          end
                      end
                  end
                  plot([0 1],[0 1],'k')
                  subplot(6,4,22);hold on;
                  plot(alls2([4 6 10 12],1),alls2([4 6 10 12],2),'ro','Markersize',5)
                  plot(alls2([4 10 12],1),alls2([4 10 12],3),'bo','Markersize',5)  
                  plot([0 10],[0 10],'k')

                inds=[4 6 12 10];
                  for k=1:4
                      i=inds(k);
                      subplot(6,4,k*4-3);hold on;
                      if Exp(i).hitabove==1
                          if ~Exp(i).postlesion
                              for m=1:length(probstay(i,1).rn)
                                  if min([numstay(i,1).rn(m) numstay(i,2).rn(m)])>10
                                      if m>Exp(i).MIN-2 % if it's hit 
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ro','Markersize',5)
                                      else
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ko','Markersize',5) 
                                      end

                                      if i~=6
                                      plot(probstay(i,1).rn(m),probstay(i,3).rn(m),'bo','Markersize',5) 
                                      end
                                  end
                              end
                          end
                      end
                        plot([0 1],[0 1],'k')
                  end

                 inds=[4 6 12 10];
                  for k=1:4
                      i=inds(k);
                      subplot(6,4,k*4-2);hold on; 
                      plot(alls2(i,1),alls2(i,2),'ro','Markersize',5)
                      if i~=6
                      plot(alls2(i,1),alls2(i,3),'bo','Markersize',5) 
                      end
                      plot([0 10],[0 10],'k')
                  end


                  subplot(6,4,19);hold on;
                  for i=indPost
                      if Exp(i).hitabove==1
                          if Exp(i).postlesion
                              for m=1:length(probstay(i,1).rn)
                                  if min([numstay(i,1).rn(m) numstay(i,2).rn(m)])>10 % if enough data
                                      if m>Exp(i).MIN-2 % if it's hit 
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ro','Markersize',5)
                                      else
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ko','Markersize',5) 
                                      end
                                  end
                              end
                          end
                      end
                  end
                 plot([0 1],[0 1],'k')

                   subplot(6,4,20);hold on;
                  plot(alls2(indPost,1),alls2(indPost,2),'ro','Markersize',5)
                  plot([0 10],[0 10],'k')


                   subplot(6,4,23);hold on;
                  for i=indPost
                      if Exp(i).hitabove==1
                          if Exp(i).postlesion
                              for m=1:length(probstay(i,1).rn)
                                  if min([numstay(i,1).rn(m) numstay(i,2).rn(m)])>10 % if enough data
                                      if m>Exp(i).MIN-2 % if it's hit 
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ro','Markersize',5)
                                      else
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ko','Markersize',5) 
                                      end
                                      plot(probstay(i,1).rn(m),probstay(i,3).rn(m),'bo','Markersize',5) 
                                  end
                              end
                          end
                      end
                  end
                 plot([0 1],[0 1],'k')

                    subplot(6,4,24);hold on;
                  plot(alls2(indPost,1),alls2(indPost,2),'ro','Markersize',5)
                  plot(alls2(indPost,1),alls2(indPost,3),'bo','Markersize',5)  
                  plot([0 10],[0 10],'k')

                inds=indPost;
                  for k=1:4
                      i=inds(k);
                      subplot(6,4,k*4-1);hold on;
                      if Exp(i).hitabove==1
                          if Exp(i).postlesion
                              for m=1:length(probstay(i,1).rn)
                                  if min([numstay(i,1).rn(m) numstay(i,2).rn(m)])>10 % if enough data
                                      if m>Exp(i).MIN-2 % if it's hit 
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ro','Markersize',5)
                                      else
                                        plot(probstay(i,1).rn(m),probstay(i,2).rn(m),'ko','Markersize',5) 
                                      end
                                     plot(probstay(i,1).rn(m),probstay(i,3).rn(m),'bo','Markersize',5) 
                                  end
                              end
                          end
                      end
                       plot([0 1],[0 1],'k')
                  end

                 inds=indPost;
                  for k=1:4
                      i=inds(k);
                      subplot(6,4,k*4);hold on; 
                      plot(alls2(i,1),alls2(i,2),'ro','Markersize',5)
                      if i~=6
                      plot(alls2(i,1),alls2(i,3),'bo','Markersize',5) 
                      end
                      plot([0 10],[0 10],'k')
                  end



