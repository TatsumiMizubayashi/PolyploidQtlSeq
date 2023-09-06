namespace PolyploidQtlSeqCore.Share
{
    /// <summary>
    /// 使用するスレッド数
    /// </summary>
    internal class ThreadNumber
    {
        /// <summary>
        /// 使用するスレッド数の最小値
        /// </summary>
        private const int MINIMUM = 1;

        /// <summary>
        /// 使用するスレッド数の最大値
        /// </summary>
        private const int MAXIMUM = 50;

        /// <summary>
        /// 使用するスレッド数を作成する。
        /// </summary>
        /// <param name="number">スレッド数</param>
        public ThreadNumber(int number)
        {
            if (number < MINIMUM || number > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(number));

            Value = number;
        }

        /// <summary>
        /// 使用するスレッド数を取得する。
        /// </summary>
        internal int Value { get; }
    }
}
