namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 3'末端トリム時のウインドウサイズ
    /// </summary>
    internal class CutTailWindowSize
    {
        /// <summary>
        /// 3'末端トリム時のウインドウサイズの最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 3'末端トリム時のウインドウサイズの最大値
        /// </summary>
        private const int MAXIMUM = 100;

        /// <summary>
        /// 3'末端トリム時のウインドウサイズを作成する。
        /// </summary>
        /// <param name="windowSize">ウインドウサイズ</param>
        public CutTailWindowSize(int windowSize)
        {
            if (windowSize < MINIMUM || windowSize > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(windowSize));

            Value = windowSize;
        }

        /// <summary>
        /// 3'末端トリム時のウインドウサイズを取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns></returns>
        internal string ToFastpArg()
        {
            return $"--cut_tail_window_size {Value}";
        }

    }
}
