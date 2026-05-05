
#V1
df2 <- data.frame(frame = V1$Frame, dist = V1$list.proj.V1) %>%
    mutate(line_color = ifelse(V11$Frame > 517, "cyan","pink"))

frame = V1$Frame
dist = V1$list.proj.V1

ggplot(df2, mapping = aes(x = frame, y = dist, color = line_color)) +
    geom_segment(aes(x = frame,  xend = lead(frame), y = dist, yend = lead(dist))) + 
    scale_color_identity()+
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = 'white', colour = 'white'),
    ) +
    xlab("Frame") +
    ylab("Distance") +
    ylim(-40,40) +
    ggtitle("V1")

#V11 <- as.data.frame(list.proj.V11)
#V11$Frame <- seq.int(nrow(V11))

df2 <- data.frame(frame = V11$Frame, dist = V11$list.proj.V11) %>%
    mutate(line_color = ifelse(V11$Frame > 0, "cyan","pink"))

frame = V11$Frame
dist = V11$list.proj.V11
           
ggplot(df2, mapping = aes(x = frame, y = dist, color = line_color)) +
    geom_segment(aes(x = frame,  xend = lead(frame), y = dist, yend = lead(dist))) + 
    scale_color_identity()+
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = 'white', colour = 'white'),
          ) +
    xlab("Frame") +
    ylab("Distance") +
    ylim(-40,40) +
    ggtitle("V11")


V21 <- as.data.frame(list.proj.V21)
V21$Frame <- seq.int(nrow(V21))


df2 <- data.frame(frame = V21$Frame, dist = V21$list.proj.V21) %>%
    mutate(line_color = ifelse(V21$Frame > 155, "cyan","pink"))

frame = V21$Frame
dist = V21$list.proj.V21

ggplot(df2, mapping = aes(x = frame, y = dist, color = line_color)) +
    geom_segment(aes(x = frame,  xend = lead(frame), y = dist, yend = lead(dist))) + 
    scale_color_identity()+
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = 'white', colour = 'white'),
    ) +
    xlab("Frame") +
    ylab("Distance") +
    ylim(-40,40) +
    ggtitle("V21")
