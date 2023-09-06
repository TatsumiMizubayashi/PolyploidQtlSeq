namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// Window移動量
    /// </summary>
    internal class StepSize
    {
        /// <summary>
        /// step sizeの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// step sizeの最大値
        /// </summary>
        private const int MAXIMUM = 100000;

        /// <summary>
        /// Stepサイズを作成する。
        /// </summary>
        /// <param name="kbpSize">step size(kbp)</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public StepSize(int kbpSize)
        {
            if (kbpSize < MINIMUM || kbpSize > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(kbpSize));

            KbpValue = kbpSize;
            BpValue = KbpValue * 1000;
        }

        /// <summary>
        /// StepSize(kbp)を取得する。
        /// </summary>
        public int KbpValue { get; }

        /// <summary>
        /// StepSize(bp)を取得する。
        /// </summary>
        public int BpValue { get; }

    }
}
