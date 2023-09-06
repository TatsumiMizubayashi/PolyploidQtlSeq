namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// Windowサイズ
    /// </summary>
    internal class WindowSize
    {
        /// <summary>
        /// window sizeの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// window sizeの最大値
        /// </summary>
        private const int MAXIMUM = 100000;

        /// <summary>
        /// WindowSizeを作成する。
        /// </summary>
        /// <param name="kbpSize">サイズ(kbp)</param>
        public WindowSize(int kbpSize)
        {
            if (kbpSize < MINIMUM || kbpSize > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(kbpSize));

            KbpValue = kbpSize;
            BpValue = KbpValue * 1000;
        }

        /// <summary>
        /// WindowSize(kbp)を取得する。
        /// </summary>
        public int KbpValue { get; }

        /// <summary>
        /// WindowSize(bp)を取得する。
        /// </summary>
        public int BpValue { get; }
    }
}
