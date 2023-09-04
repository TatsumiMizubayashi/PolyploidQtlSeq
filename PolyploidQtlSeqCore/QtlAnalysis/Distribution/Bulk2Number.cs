﻿using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// Bulk2個体数
    /// </summary>
    public class Bulk2Number
    {
        /// <summary>
        /// 個体数の最小値
        /// </summary>
        private const int MINIMUM = 2;

        /// <summary>
        /// 個体数の最大値
        /// </summary>
        private const int MAXIMUM = 1000;

        /// <summary>
        /// 個体数の規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 20;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "n2";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "NBulk2";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Number of Individuals in bulk2.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -n2 option must be an integer greater than or equal to 2 and less than or equal to 1000.";

        /// <summary>
        /// Bulk2個体数を作成する。
        /// </summary>
        /// <param name="number">個体数</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public Bulk2Number(int number, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, number, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// Bulk2個体数を作成する。
        /// </summary>
        /// <param name="number">個体数</param>
        public Bulk2Number(int number)
        {
            if (number < MINIMUM || number > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(number));

            Value = number;
        }

        /// <summary>
        /// Bulk2個体数を取得する。
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
