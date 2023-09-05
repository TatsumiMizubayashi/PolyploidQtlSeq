using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEff MaxHeap
    /// </summary>
    public class SnpEffMaxHeap
    {
        /// <summary>
        /// max heapの最小値
        /// </summary>
        private const int MINIMUM = 2;

        /// <summary>
        /// max heapの最大値
        /// </summary>
        private const int MAXIMUM = 100;

        /// <summary>
        /// max heapの規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 6;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "sm";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "snpEffMaxHeap";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "SnpEff maximum heap size (GB).";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -sm option must be an integer greater than or equal to 2 and less than or equal to 100.";

        /// <summary>
        /// SnpEff MaxHeapを作成する。
        /// </summary>
        /// <param name="maxHeap">最大ヒープサイズ(GB)</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public SnpEffMaxHeap(int maxHeap, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, maxHeap, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// SnpEff MaxHeapを作成する。
        /// </summary>
        /// <param name="maxHeap">最大ヒープサイズ(GB)</param>
        public SnpEffMaxHeap(int maxHeap)
        {
            if (maxHeap < MINIMUM || maxHeap > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(maxHeap));

            Value = maxHeap;
        }



        /// <summary>
        /// 最大ヒープサイズ(GB)を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }

    }
}
